//! Utilities related to alignment of sequences.

use anyhow::bail;
use noodles::sam::record::{cigar::op::Kind, sequence::Base, Cigar, MappingQuality};

use super::cigar::consumes_reference;
use super::cigar::consumes_sequence;

/// Filter an alignment record by its mapping quality. `true` means "filter the record" and `false` means "do not filter the record".
pub fn filter_by_mapq(
    record: &noodles::sam::alignment::Record,
    min_mapq: Option<MappingQuality>,
) -> bool {
    match min_mapq {
        Some(min_mapq) => match record.mapping_quality() {
            Some(mapq) => mapq.get() < min_mapq.get(),
            None => true, // Filter if no MAPQ is present
        },
        None => false, // Do not filter if no min MAPQ is specified
    }
}

/// Turns a condensed Cigar representation into a flattened representation. For
/// example, 10M will become a Vec of length 10 comprised completely of
/// Kind::MATCH. This utility is useful for generating a representation of a
/// Cigar string that we can step through alongside a sequence.
pub fn flatten(cigar: Cigar) -> Vec<Kind> {
    let mut result = Vec::new();

    for op in cigar.iter() {
        let len = op.len();
        let kind = op.kind();
        for _ in 1..=len {
            result.push(kind)
        }
    }

    result
}

/// Utility struct for stepping through a reference sequence, and record
/// sequence, and a Cigar string in unison.
pub struct ReferenceRecordStepThrough<'a> {
    reference_seq: &'a [Base],
    record_seq: &'a [Base],
    cigar: Vec<Kind>,
}

impl<'a> ReferenceRecordStepThrough<'a> {
    /// Creates a new [`ReferenceRecordStepThrough`].
    pub fn new(reference_seq: &'a [Base], record_seq: &'a [Base], cigar: Cigar) -> Self {
        Self {
            reference_seq,
            record_seq,
            cigar: flatten(cigar),
        }
    }

    /// Steps through a [`ReferenceRecordStepThrough`] and calls the [`Fn`] provided on
    /// the reference sequence, the record sequence (sometimes called "query" sequence),
    /// and the CIGAR string in tandem.
    pub fn stepthrough<F>(&self, mut f: F) -> anyhow::Result<()>
    where
        F: FnMut(Kind, Option<Base>, usize, Option<Base>, usize) -> anyhow::Result<()>,
    {
        let mut record_ptr = 0;
        let mut reference_ptr = 0;

        for kind in self.cigar.iter().copied() {
            let consumes_reference = consumes_reference(kind);
            let consumes_sequence = consumes_sequence(kind);

            let reference_base = match consumes_reference {
                true => {
                    let base = self.reference_seq.get(reference_ptr);
                    if let Some(base) = base {
                        Some(*base)
                    } else {
                        bail!(
                            "malformed record: record specifies that we should be able to \
                        consume a reference base, but no such base was found"
                        )
                    }
                }
                false => None,
            };

            let record_base = match consumes_sequence {
                true => {
                    let base = self.record_seq.get(record_ptr);
                    if let Some(base) = base {
                        Some(*base)
                    } else {
                        bail!(
                            "malformed record: record specifies that we should be able to \
                        consume a record base, but no such base was found"
                        )
                    }
                }
                false => None,
            };

            f(kind, reference_base, reference_ptr, record_base, record_ptr)?;

            if consumes_reference {
                reference_ptr += 1;
            }

            if consumes_sequence {
                record_ptr += 1;
            }
        }

        if self.reference_seq.len() != reference_ptr {
            bail!("reference sequence was not fully consumed");
        } else if self.record_seq.len() != record_ptr {
            bail!("record sequence was not fully consumed");
        }

        Ok(())
    }

    /// Calculates the number of edits in the [`ReferenceRecordStepThrough`] by
    /// stepping through the genome and counting up all of the mismatched M's.
    /// Errors can occur if the reference or the sequence are not all the way
    /// consumed.
    pub fn edits(&self) -> anyhow::Result<usize> {
        let mut edits = 0;

        self.stepthrough(|cigar, reference, _, record, _| {
            if cigar == Kind::Match && reference != record {
                edits += 1;
            }

            Ok(())
        })?;

        Ok(edits)
    }
}

#[cfg(test)]
mod tests {
    use noodles::sam::record::{Cigar, MappingQuality, Sequence};

    use super::*;

    #[test]
    pub fn it_filters_by_mapq() -> anyhow::Result<()> {
        let mut record = noodles::sam::alignment::Record::default();
        assert!(filter_by_mapq(
            &record,
            Some(MappingQuality::new(0).unwrap())
        )); // Get filtered because MAPQ is missing
        assert!(!filter_by_mapq(&record, None)); // Do not get filtered because filter is disabled

        record
            .mapping_quality_mut()
            .replace(MappingQuality::new(10).unwrap());
        assert!(!filter_by_mapq(
            &record,
            Some(MappingQuality::new(0).unwrap())
        )); // Do not get filtered because MAPQ is present
        assert!(!filter_by_mapq(
            &record,
            Some(MappingQuality::new(1).unwrap())
        )); // Do not get filtered because MAPQ is greater than 1
        assert!(filter_by_mapq(
            &record,
            Some(MappingQuality::new(11).unwrap())
        )); // Do get filtered because MAPQ is less than 11
        Ok(())
    }

    #[test]
    pub fn it_correctly_returns_zero_edits_when_sequences_are_identical() -> anyhow::Result<()> {
        let reference = "ACTG".parse::<Sequence>()?;
        let record = "ACTG".parse::<Sequence>()?;
        let cigar = "4M".parse::<Cigar>()?;

        let rrs = ReferenceRecordStepThrough::new(reference.as_ref(), record.as_ref(), cigar);
        assert_eq!(rrs.edits().unwrap(), 0);
        Ok(())
    }

    #[test]
    pub fn it_correctly_returns_one_edit_when_sequences_have_one_edit() -> anyhow::Result<()> {
        let reference = "AATG".parse::<Sequence>()?;
        let record = "ACTG".parse::<Sequence>()?;
        let cigar = "4M".parse::<Cigar>()?;

        let rrs = ReferenceRecordStepThrough::new(reference.as_ref(), record.as_ref(), cigar);
        assert_eq!(rrs.edits().unwrap(), 1);
        Ok(())
    }

    #[test]
    pub fn it_correctly_handles_softclips() -> anyhow::Result<()> {
        let reference = "ACTG".parse::<Sequence>()?;
        let record = "ACTGACTG".parse::<Sequence>()?;
        let cigar = "4M4S".parse::<Cigar>()?;

        let rrs = ReferenceRecordStepThrough::new(reference.as_ref(), record.as_ref(), cigar);
        assert_eq!(rrs.edits().unwrap(), 0);
        Ok(())
    }

    #[test]
    pub fn it_catches_a_malformed_record_with_too_few_record_bases() -> anyhow::Result<()> {
        let reference = "ACTG".parse::<Sequence>()?;
        let record = "ACTGACTG".parse::<Sequence>()?;
        let cigar = "4M5S".parse::<Cigar>()?;

        let rrs = ReferenceRecordStepThrough::new(reference.as_ref(), record.as_ref(), cigar);
        let result = rrs.edits();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert_eq!(
            err.to_string(),
            "malformed record: record specifies that we should \
             be able to consume a record base, but no such base was found"
        );

        Ok(())
    }

    #[test]
    pub fn it_catches_a_malformed_record_with_too_few_reference_bases() -> anyhow::Result<()> {
        let reference = "ACTG".parse::<Sequence>()?;
        let record = "ACT".parse::<Sequence>()?;
        let cigar = "3M2D".parse::<Cigar>()?;

        let rrs = ReferenceRecordStepThrough::new(reference.as_ref(), record.as_ref(), cigar);
        let result = rrs.edits();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert_eq!(
            err.to_string(),
            "malformed record: record specifies that we should \
             be able to consume a reference base, but no such base was found"
        );
        Ok(())
    }
}
