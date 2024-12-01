//! Module holding the logic for annotating junctions.

use anyhow::{bail, Ok};
use noodles::core::Position;
use noodles::sam::alignment::Record;
use noodles::sam::record::cigar::op::Kind;
use noodles::sam::record::MappingQuality;
use noodles::sam::Header;
use std::collections::{HashMap, HashSet};

use crate::derive::junction_annotation::results;
use crate::utils::alignment::filter_by_mapq;

/// Struct to hold starts and ends of exons.
pub struct ExonSets<'a> {
    /// Starts of exons, grouped by contig.
    pub starts: HashMap<&'a str, HashSet<Position>>,

    /// ends of exons, grouped by contig.
    pub ends: HashMap<&'a str, HashSet<Position>>,
}

/// Parameters defining how to annotate found junctions
pub struct JunctionAnnotationParameters {
    /// Minimum intron length to consider.
    pub min_intron_length: usize,

    /// Minimum number of reads supporting a junction to be considered.
    pub min_read_support: usize,

    /// Minumum mapping quality for a record to be considered.
    /// `None` means no filtering by MAPQ. This also allows
    /// for records _without_ a MAPQ to be counted.
    pub min_mapq: Option<MappingQuality>,

    /// Do not count supplementary alignments.
    pub no_supplementary: bool,

    /// Do count secondary alignments.
    pub count_secondary: bool,

    /// Do count duplicates.
    pub count_duplicates: bool,
}

/// Function for incrementing a junction counter by one.
fn increment_junction_counter(
    junction_counter: &mut results::JunctionCounter,
    junction: results::Junction,
) {
    junction_counter
        .entry(junction)
        .and_modify(|e| *e += 1)
        .or_insert(1);
}

/// Function for incrementing a junction map by one.
fn increment_junction_map(
    junction_map: &mut results::JunctionsMap,
    ref_name: &str,
    junction: results::Junction,
) {
    increment_junction_counter(
        junction_map.entry(ref_name.to_string()).or_default(),
        junction,
    );
}

/// Function to filter out records based on their flags.
fn filter_by_flags(record: &Record, params: &JunctionAnnotationParameters) -> bool {
    let flags = record.flags();
    if flags.is_unmapped()
        || (params.no_supplementary && flags.is_supplementary())
        || (!params.count_secondary && flags.is_secondary())
        || (!params.count_duplicates && flags.is_duplicate())
    {
        return true;
    }
    false
}

/// Function to filter out records that don't have introns.
fn filter_by_cigar(record: &Record) -> bool {
    !record
        .cigar()
        .iter()
        .any(|op| matches!(op.kind(), Kind::Skip))
}

/// Main function to annotate junctions one record at a time.
pub fn process(
    record: &Record,
    exons: &ExonSets<'_>,
    header: &Header,
    params: &JunctionAnnotationParameters,
    results: &mut results::JunctionAnnotationResults,
) -> anyhow::Result<()> {
    // (1) Parse the read name.
    let read_name = match record.read_name() {
        Some(name) => name,
        _ => bail!("Could not parse read name"),
    };

    // (2) Filter by record flags.
    if filter_by_flags(record, params) {
        results.records.filtered_by_flags += 1;
        return Ok(());
    }

    // (3) Filter by CIGAR.
    // We only care about reads with introns, so if there are no introns
    // we can skip this read.
    if filter_by_cigar(record) {
        results.records.not_spliced += 1;
        return Ok(());
    }

    // (4) Filter by MAPQ
    if filter_by_mapq(record, params.min_mapq) {
        results.records.bad_mapq += 1;
        return Ok(());
    }

    // (5) Parse the reference sequence from the record.
    let (seq_name, _) = match record.reference_sequence(header) {
        Some(seq_map_result) => seq_map_result?,
        _ => {
            bail!(
                "Could not parse reference sequence id for read: {}",
                read_name
            )
        }
    };
    let seq_name = seq_name.as_str();

    // (6) Check if there will be annotations for this reference sequence.
    let mut ref_is_annotated = true;
    if !exons.starts.contains_key(seq_name) || !exons.ends.contains_key(seq_name) {
        ref_is_annotated = false;
    }

    // (7) Calculate the start position of this read. This will
    // be used to find the position of any introns.
    let start = match record.alignment_start() {
        Some(s) => s,
        _ => bail!("Could not parse record's start position."),
    };

    // (8) Find introns
    let mut cur_pos = start;
    for op in record.cigar().iter() {
        match op.kind() {
            // This is an intron.
            Kind::Skip => {
                // Check that `op.len() >= params.min_intron_length` later,
                // once all reads supporting short junctions have been collected
                // for better metric reporting.

                let intron_start = cur_pos;
                // Update cur_pos to the end of the intron.
                cur_pos = cur_pos.checked_add(op.len()).unwrap();
                let intron_end = cur_pos;
                let junction: results::Junction = (intron_start, intron_end);

                // If the reference sequence is not annotated, we can skip
                // the lookup of exon positions, and directly insert the
                // intron into the unannotated_reference HashMap.
                if !ref_is_annotated {
                    increment_junction_map(
                        &mut results.junction_annotations.unannotated_reference,
                        seq_name,
                        junction,
                    );
                    continue;
                }

                // The following unwraps are safe because we checked that the reference
                // sequence is annotated above.
                let exon_starts = exons.starts.get(seq_name).unwrap();
                let exon_ends = exons.ends.get(seq_name).unwrap();

                let mut intron_start_known = false;
                let mut intron_end_known = false;
                if exon_ends.contains(&intron_start) {
                    intron_start_known = true;
                }
                if exon_starts.contains(&intron_end) {
                    intron_end_known = true;
                }

                let junction_map = match (intron_start_known, intron_end_known) {
                    (true, true) => {
                        // We found both ends of the intron.
                        // This is a Known Junction.
                        &mut results.junction_annotations.known
                    }
                    (true, false) | (false, true) => {
                        // We found one end of the intron,
                        // but not the other.
                        // This is a Partial Novel Junction.
                        &mut results.junction_annotations.partial_novel
                    }
                    (false, false) => {
                        // We found neither end of the intron.
                        // This is a Complete Novel Junction.
                        &mut results.junction_annotations.complete_novel
                    }
                };
                increment_junction_map(junction_map, seq_name, junction)
            }
            // Operations that increment the reference position (beside Skip which is handled above).
            Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                cur_pos = cur_pos.checked_add(op.len()).unwrap();
            }
            // Operations that do not increment the reference position.
            _ => {}
        }
    }

    results.records.processed += 1;
    Ok(())
}

/// Function to filter out junctions that are too short or don't have enough read support.
fn filter_junction_map(
    junction_map: &mut results::JunctionsMap,
    min_intron_length: usize,
    min_read_support: usize,
    metrics: &mut results::SummaryResults,
) {
    junction_map.retain(|_, v| {
        v.retain(|(start, end), count| {
            let mut keep = true;
            if end.get() - start.get() < min_intron_length {
                metrics.intron_too_short += 1;
                keep = false;
            }
            if *count < min_read_support {
                metrics.junctions_with_not_enough_read_support += 1;
                keep = false;
            }
            if !keep {
                metrics.total_rejected_junctions += 1;
            }
            keep
        });
        !v.is_empty()
    });
}

/// Function to tally up the junctions and their read support.
fn tally_junctions_and_support(junction_map: &results::JunctionsMap) -> (usize, usize) {
    let junctions = junction_map.values().map(|v| v.len()).sum();
    let support = junction_map
        .values()
        .map(|v| v.values().sum::<usize>())
        .sum();
    (junctions, support)
}

/// Main function to summarize the results of the junction-annotation subcommand.
pub fn summarize(
    results: &mut results::JunctionAnnotationResults,
    params: &JunctionAnnotationParameters,
) {
    // Filter out junctions that are too short or don't have enough read support.
    filter_junction_map(
        &mut results.junction_annotations.known,
        params.min_intron_length,
        params.min_read_support,
        &mut results.summary,
    );
    filter_junction_map(
        &mut results.junction_annotations.partial_novel,
        params.min_intron_length,
        params.min_read_support,
        &mut results.summary,
    );
    filter_junction_map(
        &mut results.junction_annotations.complete_novel,
        params.min_intron_length,
        params.min_read_support,
        &mut results.summary,
    );
    filter_junction_map(
        &mut results.junction_annotations.unannotated_reference,
        params.min_intron_length,
        params.min_read_support,
        &mut results.summary,
    );

    // Tally up observed junctions and spliced reads.
    (
        results.summary.known_junctions,
        results.summary.known_junctions_read_support,
    ) = tally_junctions_and_support(&results.junction_annotations.known);
    (
        results.summary.partial_novel_junctions,
        results.summary.partial_novel_junctions_read_support,
    ) = tally_junctions_and_support(&results.junction_annotations.partial_novel);
    (
        results.summary.complete_novel_junctions,
        results.summary.complete_novel_junctions_read_support,
    ) = tally_junctions_and_support(&results.junction_annotations.complete_novel);
    (
        results.summary.unannotated_reference_junctions,
        results.summary.unannotated_reference_junctions_read_support,
    ) = tally_junctions_and_support(&results.junction_annotations.unannotated_reference);

    // Tally up total junctions.
    results.summary.total_junctions = results.summary.known_junctions
        + results.summary.partial_novel_junctions
        + results.summary.complete_novel_junctions
        + results.summary.unannotated_reference_junctions;
    // Tally up total read support.
    results.summary.total_junctions_read_support = results.summary.known_junctions_read_support
        + results.summary.partial_novel_junctions_read_support
        + results.summary.complete_novel_junctions_read_support
        + results.summary.unannotated_reference_junctions_read_support;

    // Calculate percentages.
    let total_junctions = results.summary.total_junctions as f64
        - results.summary.unannotated_reference_junctions as f64; // exclude unannotated junctions from percentages
    results.summary.known_junctions_percent =
        results.summary.known_junctions as f64 / total_junctions * 100.0;
    results.summary.partial_novel_junctions_percent =
        results.summary.partial_novel_junctions as f64 / total_junctions * 100.0;
    results.summary.complete_novel_junctions_percent =
        results.summary.complete_novel_junctions as f64 / total_junctions * 100.0;

    // Calculate average read support.
    // Total
    results.summary.average_junction_read_support = results.summary.total_junctions_read_support
        as f64
        / results.summary.total_junctions as f64;
    // Known
    results.summary.average_known_junction_read_support =
        results.summary.known_junctions_read_support as f64
            / results.summary.known_junctions as f64;
    // Partial Novel
    results.summary.average_partial_novel_junction_read_support =
        results.summary.partial_novel_junctions_read_support as f64
            / results.summary.partial_novel_junctions as f64;
    // Complete Novel
    results.summary.average_complete_novel_junction_read_support =
        results.summary.complete_novel_junctions_read_support as f64
            / results.summary.complete_novel_junctions as f64;
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::core::Position;
    use noodles::sam::header::record::value::map;
    use noodles::sam::header::record::value::map::header::Version;
    use noodles::sam::header::record::value::map::{Map, ReferenceSequence};
    use noodles::sam::record::MappingQuality;
    use noodles::sam::record::ReadName;
    use std::num::NonZeroUsize;

    fn create_test_exons() -> ExonSets<'static> {
        let exon_starts: HashMap<&str, HashSet<Position>> = HashMap::from([(
            "sq1",
            HashSet::from([
                Position::new(1).unwrap(),
                Position::new(11).unwrap(),
                Position::new(21).unwrap(),
                Position::new(31).unwrap(),
                Position::new(41).unwrap(),
                Position::new(51).unwrap(),
                Position::new(61).unwrap(),
                Position::new(71).unwrap(),
            ]),
        )]);
        let exon_ends: HashMap<&str, HashSet<Position>> = exon_starts
            .iter()
            .map(|(k, v)| (*k, v.iter().map(|e| e.checked_add(10).unwrap()).collect()))
            .collect::<HashMap<&str, HashSet<Position>>>();
        let exons: ExonSets<'_> = ExonSets {
            starts: exon_starts,
            ends: exon_ends,
        };
        exons
    }

    fn create_test_header() -> Header {
        Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence(
                "sq1".parse().unwrap(),
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(800).unwrap()),
            )
            .build()
    }

    #[test]
    fn test_filter_by_flags() {
        // Setup
        let mut record = Record::default();
        let params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
            no_supplementary: false,
            count_secondary: false,
            count_duplicates: false,
        };

        // Test that records are filtered out correctly
        record.flags_mut().set(0x4.into(), true);
        assert!(filter_by_flags(&record, &params));
        record.flags_mut().set(0x4.into(), false);
        record.flags_mut().set(0x800.into(), true);
        assert!(!filter_by_flags(&record, &params));
        record.flags_mut().set(0x800.into(), false);
        record.flags_mut().set(0x100.into(), true);
        assert!(filter_by_flags(&record, &params));
        record.flags_mut().set(0x100.into(), false);
        record.flags_mut().set(0x400.into(), true);
        assert!(filter_by_flags(&record, &params));
        record.flags_mut().set(0x400.into(), false);
        assert!(!filter_by_flags(&record, &params));
    }

    #[test]
    fn test_filter_by_cigar() {
        // Setup
        let mut record = Record::default();

        // Test that records are filtered out correctly
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        assert!(!filter_by_cigar(&record));
        *record.cigar_mut() = "10M".parse().unwrap();
        assert!(filter_by_cigar(&record));
    }

    #[test]
    fn test_filter_junction_map() {
        // Setup
        let mut junction_map = results::JunctionsMap::default();
        junction_map.insert(
            "sq1".to_string(),
            HashMap::from([
                ((Position::new(1).unwrap(), Position::new(10).unwrap()), 3),
                ((Position::new(1).unwrap(), Position::new(11).unwrap()), 1),
                ((Position::new(1).unwrap(), Position::new(5).unwrap()), 1),
            ]),
        );
        junction_map.insert(
            "sq2".to_string(),
            HashMap::from([((Position::new(1).unwrap(), Position::new(11).unwrap()), 2)]),
        );
        let min_intron_length = 10;
        let min_read_support = 2;
        let mut metrics = results::SummaryResults::default();

        // Test that junctions are filtered out correctly
        filter_junction_map(
            &mut junction_map,
            min_intron_length,
            min_read_support,
            &mut metrics,
        );
        assert_eq!(junction_map.len(), 1);
        assert_eq!(junction_map.get("sq1"), None);
        assert_eq!(junction_map.get("sq2").unwrap().len(), 1);
        assert_eq!(metrics.intron_too_short, 2);
        assert_eq!(metrics.junctions_with_not_enough_read_support, 2);
        assert_eq!(metrics.total_rejected_junctions, 3);
    }

    #[test]
    fn test_tally_junctions_and_support() {
        // Setup
        let mut junction_map = results::JunctionsMap::default();
        junction_map.insert(
            "sq1".to_string(),
            HashMap::from([
                ((Position::new(1).unwrap(), Position::new(11).unwrap()), 1),
                ((Position::new(1).unwrap(), Position::new(5).unwrap()), 1),
            ]),
        );
        junction_map.insert(
            "sq2".to_string(),
            HashMap::from([((Position::new(1).unwrap(), Position::new(11).unwrap()), 2)]),
        );

        // Test that junctions are tallied correctly
        let (juncs, support) = tally_junctions_and_support(&junction_map);
        assert_eq!(juncs, 3);
        assert_eq!(support, 4);
    }

    #[test]
    fn test_process_known_junction() {
        // Setup
        let mut results = results::JunctionAnnotationResults::default();
        let params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
            no_supplementary: false,
            count_secondary: false,
            count_duplicates: false,
        };
        let exons = create_test_exons();
        let header = create_test_header();

        // Test known junction
        let mut record = Record::default();
        let r1_name: ReadName = "known1".parse().unwrap();
        *record.read_name_mut() = Some(r1_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 1);
        assert_eq!(results.records.filtered_by_flags, 0);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);
        assert_eq!(results.junction_annotations.known.len(), 1);
        assert_eq!(
            results.junction_annotations.known.get("sq1").unwrap().len(),
            1
        );
        assert_eq!(
            results
                .junction_annotations
                .known
                .get("sq1")
                .unwrap()
                .get(&(Position::new(11).unwrap(), Position::new(21).unwrap()))
                .unwrap(),
            &1
        );
    }

    #[test]
    fn test_process_partial_novel_junction() {
        // Setup
        let mut results = results::JunctionAnnotationResults::default();
        let params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
            no_supplementary: false,
            count_secondary: false,
            count_duplicates: false,
        };
        let exons = create_test_exons();
        let header = create_test_header();

        // Test partial novel junction
        let mut record = Record::default();
        let r1_name: ReadName = "partial1".parse().unwrap();
        *record.read_name_mut() = Some(r1_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M12N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 1);
        assert_eq!(results.records.filtered_by_flags, 0);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);
        assert_eq!(results.junction_annotations.partial_novel.len(), 1);
        assert_eq!(
            results
                .junction_annotations
                .partial_novel
                .get("sq1")
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            results
                .junction_annotations
                .partial_novel
                .get("sq1")
                .unwrap()
                .get(&(Position::new(11).unwrap(), Position::new(23).unwrap()))
                .unwrap(),
            &1
        );
    }

    #[test]
    fn test_process_complete_novel_junction() {
        // Setup
        let mut results = results::JunctionAnnotationResults::default();
        let params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
            no_supplementary: false,
            count_secondary: false,
            count_duplicates: false,
        };
        let exons = create_test_exons();
        let header = create_test_header();

        // Test complete novel junction
        let mut record = Record::default();
        let r1_name: ReadName = "complete1".parse().unwrap();
        *record.read_name_mut() = Some(r1_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "85M14N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 1);
        assert_eq!(results.records.filtered_by_flags, 0);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);
        assert_eq!(results.junction_annotations.complete_novel.len(), 1);
        assert_eq!(
            results
                .junction_annotations
                .complete_novel
                .get("sq1")
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            results
                .junction_annotations
                .complete_novel
                .get("sq1")
                .unwrap()
                .get(&(Position::new(86).unwrap(), Position::new(100).unwrap()))
                .unwrap(),
            &1
        );
    }

    #[test]
    fn test_process_ignores_unmapped() {
        // Setup
        let mut results = results::JunctionAnnotationResults::default();
        let params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
            no_supplementary: false,
            count_secondary: false,
            count_duplicates: false,
        };
        let exons = create_test_exons();
        let header = create_test_header();

        // Test that unmapped gets ignored
        let mut record = Record::default();
        let r1_name: ReadName = "unmapped".parse().unwrap();
        *record.read_name_mut() = Some(r1_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), true);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 0);
        assert_eq!(results.records.filtered_by_flags, 1);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);
    }

    #[test]
    fn test_process_supplementary_toggle() {
        // Setup
        let mut results = results::JunctionAnnotationResults::default();
        let mut params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
            no_supplementary: true,
            count_secondary: false,
            count_duplicates: false,
        };
        let exons = create_test_exons();
        let header = create_test_header();

        // Test that supplementary gets ignored
        let mut record = Record::default();
        let r1_name: ReadName = "supplementary1".parse().unwrap();
        *record.read_name_mut() = Some(r1_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        record.flags_mut().set(0x800.into(), true);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 0);
        assert_eq!(results.records.filtered_by_flags, 1);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);

        // Test that supplementary gets processed
        params.no_supplementary = false;

        let mut record = Record::default();
        let r2_name = "supplementary2".parse().unwrap();
        *record.read_name_mut() = Some(r2_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        record.flags_mut().set(0x800.into(), false);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 1);
        assert_eq!(results.records.filtered_by_flags, 1);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);
    }

    #[test]
    fn test_process_secondary_toggle() {
        // Setup
        let mut results = results::JunctionAnnotationResults::default();
        let mut params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
            no_supplementary: false,
            count_secondary: true,
            count_duplicates: false,
        };
        let exons = create_test_exons();
        let header = create_test_header();

        // Test that secondary gets processed
        let mut record = Record::default();
        let r1_name: ReadName = "secondary1".parse().unwrap();
        *record.read_name_mut() = Some(r1_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        record.flags_mut().set(0x100.into(), true);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 1);
        assert_eq!(results.records.filtered_by_flags, 0);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);

        // Test that secondary gets ignored
        params.count_secondary = false;

        let mut record = Record::default();
        let r2_name = "secondary2".parse().unwrap();
        *record.read_name_mut() = Some(r2_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        record.flags_mut().set(0x100.into(), true);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 1);
        assert_eq!(results.records.filtered_by_flags, 1);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);
    }

    #[test]
    fn test_process_mapq_toggle() {
        // Setup
        let mut results = results::JunctionAnnotationResults::default();
        let mut params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
            no_supplementary: false,
            count_secondary: false,
            count_duplicates: false,
        };
        let exons = create_test_exons();
        let header = create_test_header();

        // Test that mapq gets processed
        let mut record = Record::default();
        let r1_name: ReadName = "mapq1".parse().unwrap();
        *record.read_name_mut() = Some(r1_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 1);
        assert_eq!(results.records.filtered_by_flags, 0);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);

        // Test that mapq gets ignored
        params.min_mapq = Some(MappingQuality::new(61).unwrap());

        let mut record = Record::default();
        let r2_name = "mapq2".parse().unwrap();
        *record.read_name_mut() = Some(r2_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 1);
        assert_eq!(results.records.filtered_by_flags, 0);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 1);
    }

    #[test]
    fn test_process_intron_too_short() {
        // Setup
        let mut results = results::JunctionAnnotationResults::default();
        let params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
            no_supplementary: false,
            count_secondary: false,
            count_duplicates: false,
        };
        let exons = create_test_exons();
        let header = create_test_header();

        // Test that intron too short gets processed
        let mut record = Record::default();
        let r1_name: ReadName = "short1".parse().unwrap();
        *record.read_name_mut() = Some(r1_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M5N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 1); // processed at first, gets filtered later
        assert_eq!(results.records.filtered_by_flags, 0);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);
    }

    #[test]
    fn test_process_multiple_junctions() {
        // Setup
        let mut results = results::JunctionAnnotationResults::default();
        let params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
            no_supplementary: false,
            count_secondary: false,
            count_duplicates: false,
        };
        let exons = create_test_exons();
        let header = create_test_header();

        // Test that multiple junctions are processed
        let mut record = Record::default();
        let r1_name: ReadName = "long_read".parse().unwrap();
        *record.read_name_mut() = Some(r1_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M10N10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 1);
        assert_eq!(results.records.filtered_by_flags, 0);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);
        assert_eq!(results.junction_annotations.known.len(), 1);
        assert_eq!(
            results.junction_annotations.known.get("sq1").unwrap().len(),
            3
        );
        assert_eq!(
            results
                .junction_annotations
                .known
                .get("sq1")
                .unwrap()
                .get(&(Position::new(11).unwrap(), Position::new(21).unwrap()))
                .unwrap(),
            &1
        );
        assert_eq!(
            results
                .junction_annotations
                .known
                .get("sq1")
                .unwrap()
                .get(&(Position::new(31).unwrap(), Position::new(41).unwrap()))
                .unwrap(),
            &1
        );
        assert_eq!(
            results
                .junction_annotations
                .known
                .get("sq1")
                .unwrap()
                .get(&(Position::new(51).unwrap(), Position::new(61).unwrap()))
                .unwrap(),
            &1
        );
    }

    #[test]
    fn test_process_unspliced_read() {
        // Setup
        let mut results = results::JunctionAnnotationResults::default();
        let params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
            no_supplementary: false,
            count_secondary: false,
            count_duplicates: false,
        };
        let exons = create_test_exons();
        let header = create_test_header();

        // Test that unspliced gets ignored
        let mut record = Record::default();
        let r1_name: ReadName = "unspliced".parse().unwrap();
        *record.read_name_mut() = Some(r1_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        process(&record, &exons, &header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 0);
        assert_eq!(results.records.filtered_by_flags, 0);
        assert_eq!(results.records.not_spliced, 1);
        assert_eq!(results.records.bad_mapq, 0);
    }

    #[test]
    fn test_process_unannotated_reference() {
        // Setup
        let mut results = results::JunctionAnnotationResults::default();
        let params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
            no_supplementary: false,
            count_secondary: false,
            count_duplicates: false,
        };
        let rand_header = Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence(
                "sq1_random".parse().unwrap(),
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(800).unwrap()),
            )
            .build();
        let exons = create_test_exons();

        // Test that unannotated reference gets processed
        let mut record = Record::default();
        let r1_name: ReadName = "unannotated".parse().unwrap();
        *record.read_name_mut() = Some(r1_name);
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(1);
        *record.cigar_mut() = "10M10N10M".parse().unwrap();
        *record.mapping_quality_mut() = MappingQuality::new(60);
        record.flags_mut().set(0x4.into(), false);
        process(&record, &exons, &rand_header, &params, &mut results).unwrap();
        assert_eq!(results.records.processed, 1);
        assert_eq!(results.records.filtered_by_flags, 0);
        assert_eq!(results.records.not_spliced, 0);
        assert_eq!(results.records.bad_mapq, 0);
        assert_eq!(results.junction_annotations.unannotated_reference.len(), 1);
        assert_eq!(
            results
                .junction_annotations
                .unannotated_reference
                .get("sq1_random")
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            results
                .junction_annotations
                .unannotated_reference
                .get("sq1_random")
                .unwrap()
                .get(&(Position::new(11).unwrap(), Position::new(21).unwrap()))
                .unwrap(),
            &1
        );
    }

    #[test]
    fn test_summarize() {
        // Setup
        let mut results = results::JunctionAnnotationResults::default();
        let params = JunctionAnnotationParameters {
            min_intron_length: 10,
            min_read_support: 2,
            min_mapq: Some(MappingQuality::new(30).unwrap()),
            no_supplementary: false,
            count_secondary: false,
            count_duplicates: false,
        };
        results.junction_annotations.known.insert(
            "sq1".to_string(),
            HashMap::from([
                ((Position::new(11).unwrap(), Position::new(21).unwrap()), 4),
                ((Position::new(31).unwrap(), Position::new(41).unwrap()), 1),
                ((Position::new(21).unwrap(), Position::new(41).unwrap()), 3),
            ]),
        );
        results.junction_annotations.partial_novel.insert(
            "sq1".to_string(),
            HashMap::from([
                ((Position::new(11).unwrap(), Position::new(37).unwrap()), 3),
                ((Position::new(11).unwrap(), Position::new(15).unwrap()), 2),
            ]),
        );
        results.junction_annotations.complete_novel.insert(
            "sq1".to_string(),
            HashMap::from([(
                (Position::new(103).unwrap(), Position::new(117).unwrap()),
                2,
            )]),
        );
        results.junction_annotations.unannotated_reference.insert(
            "sq1_random".to_string(),
            HashMap::from([((Position::new(1).unwrap(), Position::new(11).unwrap()), 5)]),
        );

        // Test that results are summarized correctly
        summarize(&mut results, &params);
        assert_eq!(results.summary.known_junctions, 2);
        assert_eq!(results.summary.known_junctions_read_support, 7);
        assert_eq!(results.summary.partial_novel_junctions, 1);
        assert_eq!(results.summary.partial_novel_junctions_read_support, 3);
        assert_eq!(results.summary.complete_novel_junctions, 1);
        assert_eq!(results.summary.complete_novel_junctions_read_support, 2);
        assert_eq!(results.summary.unannotated_reference_junctions, 1);
        assert_eq!(
            results.summary.unannotated_reference_junctions_read_support,
            5
        );
        assert_eq!(results.summary.total_junctions, 5);
        assert_eq!(results.summary.total_junctions_read_support, 17);
        assert_eq!(results.summary.known_junctions_percent, 50.0);
        assert_eq!(results.summary.partial_novel_junctions_percent, 25.0);
        assert_eq!(results.summary.complete_novel_junctions_percent, 25.0);
        assert_eq!(results.summary.average_junction_read_support, 3.4);
        assert_eq!(results.summary.average_known_junction_read_support, 3.5);
        assert_eq!(
            results.summary.average_partial_novel_junction_read_support,
            3.0
        );
        assert_eq!(
            results.summary.average_complete_novel_junction_read_support,
            2.0
        );
    }
}
