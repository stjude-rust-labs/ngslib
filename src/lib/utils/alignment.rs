use noodles_sam::record::{cigar::op::Kind, sequence::Base, Cigar};

use super::cigar::{consumes_reference, consumes_sequence};

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

pub struct ReferenceRecordStepThrough<'a> {
    reference_seq: &'a [u8],
    record_seq: &'a [Base],
    cigar: Vec<Kind>,
}

impl<'a> ReferenceRecordStepThrough<'a> {
    pub fn new(reference_seq: &'a [u8], record_seq: &'a [Base], cigar: Cigar) -> Self {
        Self {
            reference_seq,
            record_seq,
            cigar: flatten(cigar),
        }
    }

    pub fn edits(&self) -> usize {
        let mut edits = 0;
        let mut record_ptr = 0;
        let mut reference_ptr = 0;

        for kind in self.cigar.iter().copied() {
            if kind == Kind::Match {
                let ref_base = self.reference_seq[reference_ptr] as char;
                let record_base: char = self.record_seq[record_ptr].into();
                if ref_base != record_base {
                    edits += 1;
                }
            }

            if consumes_reference(kind) {
                reference_ptr += 1;
            }

            if consumes_sequence(kind) {
                record_ptr += 1;
            }
        }

        edits
    }
}
