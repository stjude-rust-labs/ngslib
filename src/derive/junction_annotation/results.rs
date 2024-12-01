//! Results related to the `ngs derive junction-annotation` subcommand.

use noodles::core::Position;
use serde::ser::SerializeStruct;
use serde::{Serialize, Serializer};
use std::collections::HashMap;

/// A junction is a tuple of (start, end) coordinates.
pub type Junction = (Position, Position);

/// A junction counter is a HashMap where the key is a junction and the value is the number of
/// spliced reads that support the junction.
pub type JunctionCounter = HashMap<Junction, usize>;

/// A map of junctions. The key is the reference name, and the value is a JunctionCounter.
pub type JunctionsMap = HashMap<String, JunctionCounter>;

/// Lists of annotated junctions.
#[derive(Clone, Debug, Default)]
pub struct JunctionAnnotations {
    /// Known junctions. The outer key is the referece name, and the value is another
    /// HashMap. The inner key is the (start, end) coordinates of a junction,
    /// and the value is the number of spliced reads that support the junction.
    pub known: JunctionsMap,

    /// Partially novel junctions. The outer key is the referece name, and the value is another
    /// HashMap. The inner key is the (start, end) coordinates of a junction,
    /// and the value is the number of spliced reads that support the junction.
    pub partial_novel: JunctionsMap,

    /// Complete novel junctions. The outer key is the referece name, and the value is another
    /// HashMap. The inner key is the (start, end) coordinates of a junction,
    /// and the value is the number of spliced reads that support the junction.
    pub complete_novel: JunctionsMap,

    /// Junctions on reference sequences for which junction annotations were not found.
    /// The outer key is the referece name, and the value is another
    /// HashMap. The inner key is the (start, end) coordinates of a junction,
    /// and the value is the number of spliced reads that support the junction.
    pub unannotated_reference: JunctionsMap,
}

// TODO should contigs be sorted?
impl Serialize for JunctionAnnotations {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut known = Vec::new();
        for (ref_name, junctions) in &self.known {
            let mut junctions_vec = Vec::new();
            for ((start, end), count) in junctions {
                junctions_vec.push((start.get(), end.get(), count));
            }
            known.push((ref_name.clone(), junctions_vec));
        }

        let mut partial_novel = Vec::new();
        for (ref_name, junctions) in &self.partial_novel {
            let mut junctions_vec = Vec::new();
            for ((start, end), count) in junctions {
                junctions_vec.push((start.get(), end.get(), count));
            }
            partial_novel.push((ref_name.clone(), junctions_vec));
        }

        let mut complete_novel = Vec::new();
        for (ref_name, junctions) in &self.complete_novel {
            let mut junctions_vec = Vec::new();
            for ((start, end), count) in junctions {
                junctions_vec.push((start.get(), end.get(), count));
            }
            complete_novel.push((ref_name.clone(), junctions_vec));
        }

        let mut unannotated_reference = Vec::new();
        for (ref_name, junctions) in &self.unannotated_reference {
            let mut junctions_vec = Vec::new();
            for ((start, end), count) in junctions {
                junctions_vec.push((start.get(), end.get(), count));
            }
            unannotated_reference.push((ref_name.clone(), junctions_vec));
        }

        let mut s = serializer.serialize_struct("JunctionAnnotations", 4)?;
        s.serialize_field("known", &known)?;
        s.serialize_field("partial_novel", &partial_novel)?;
        s.serialize_field("complete_novel", &complete_novel)?;
        s.serialize_field("unannotated_reference", &unannotated_reference)?;
        s.end()
    }
}

/// General record metrics that are tallied as a part of the
/// junction-annotation subcommand.
#[derive(Clone, Debug, Default, Serialize)]
pub struct RecordMetrics {
    /// The number of records that have been fully processed.
    /// This is the number of spliced records that have been considered.
    pub processed: usize,

    /// The number of records that have been ignored because of their flags.
    /// (i.e. they were unmapped, duplicates, secondary, or supplementary)
    /// The last 3 conditions can be toggled on/off with CL flags
    pub filtered_by_flags: usize,

    /// The number of records that have been ignored because they were not
    /// spliced.
    pub not_spliced: usize,

    /// The number of records with junctions that have been ignored because
    /// they failed the MAPQ filter.
    /// This could either mean the MAPQ was too low or it was missing.
    pub bad_mapq: usize,
}

/// Summary statistics for the junction-annotation subcommand.
#[derive(Clone, Default, Debug, Serialize)]
pub struct SummaryResults {
    /// The total number of junctions observed in the file.
    pub total_junctions: usize,

    /// The total number of splices observed in the file.
    /// More than one splice can be observed per read, especially
    /// with long read data, so this number is not necessarily equal
    /// to the number of spliced reads. It may be greater.
    pub total_junctions_read_support: usize,

    /// The average number of spliced reads supporting a junction.
    pub average_junction_read_support: f64,

    /// The total number of known junctions observed in the file.
    pub known_junctions: usize,

    ///The total number of partially novel junctions observed in the file.
    pub partial_novel_junctions: usize,

    /// The total number of complete novel junctions observed in the file.
    pub complete_novel_junctions: usize,

    /// The total number of junctions on reference sequences for which junction
    /// annotations were not found.
    pub unannotated_reference_junctions: usize,

    /// The number of reads supporting known junctions.
    /// If a read supports more than one known junction, it is counted more than once.
    /// A read with more one junction may also contribute to the support of
    /// partially novel or completely novel junctions.
    pub known_junctions_read_support: usize,

    /// The number of reads supporting partially novel junctions.
    /// If a read supports more than one partially novel junction, it is counted more than once.
    /// A read with more one junction may also contribute to the support of
    /// known or completely novel junctions.
    pub partial_novel_junctions_read_support: usize,

    /// The number of reads supporting completely novel junctions.
    /// If a read supports more than one completely novel junction, it is counted more than once.
    /// A read with more one junction may also contribute to the support of
    /// known or partially novel junctions.
    pub complete_novel_junctions_read_support: usize,

    /// The number of reads supporting junctions on reference sequences for which
    /// junction annotations were not found.
    /// If a read supports more than one junction, it is counted more than once.
    pub unannotated_reference_junctions_read_support: usize,

    /// The percentage of junctions that are known.
    /// This percentage excludes junctions on reference sequences for which
    /// junction annotations were not found.
    pub known_junctions_percent: f64,

    /// The percentage of junctions that are partially novel.
    /// This percentage excludes junctions on reference sequences for which
    /// junction annotations were not found.
    pub partial_novel_junctions_percent: f64,

    /// The percentage of junctions that are completely novel.
    /// This percentage excludes junctions on reference sequences for which
    /// junction annotations were not found.
    pub complete_novel_junctions_percent: f64,

    /// Average number of reads supporting known junctions.
    pub average_known_junction_read_support: f64,

    /// Average number of reads supporting partially novel junctions.
    pub average_partial_novel_junction_read_support: f64,

    /// Average number of reads supporting completely novel junctions.
    pub average_complete_novel_junction_read_support: f64,

    /// The total number of junctions that have been rejected because
    /// they failed the --min-read-support or the --min-intron-length filter.
    /// A junction can be rejected for both reasons, so this
    /// number may not be equal to the sum of junctions_with_not_enough_read_support
    /// and intron_too_short.
    pub total_rejected_junctions: usize,

    /// The total number of junctions which were discarded due to lack of
    /// read support. This is not mutually exclusive with intron_too_short.
    pub junctions_with_not_enough_read_support: usize,

    /// The number of junctions that have been ignored because
    /// they failed the min_intron_length filter.
    /// This is not mutually exclusive with junctions_with_not_enough_read_support.
    pub intron_too_short: usize,
}

/// Main Results struct. This struct aggregates all of the minor metrics structs
/// outlined in this file so they can be kept track of as a unit.
#[derive(Clone, Default, Debug, Serialize)]
pub struct JunctionAnnotationResults {
    /// Lists of annotated junctions.
    pub junction_annotations: JunctionAnnotations,

    /// General record metrics.
    pub records: RecordMetrics,

    /// Summary statistics for the junction-annotation subcommand.
    pub summary: SummaryResults,
}
