//! Results structs for the strandedness subcommand.

use serde::Serialize;

/// General read record metrics that are tallied as a part of the
/// strandedness subcommand.
#[derive(Clone, Default, Serialize, Debug)]
pub struct ReadRecordMetrics {
    /// The number of records that have been filtered because of their flags.
    /// (i.e. they were qc_fail, duplicates, secondary, or supplementary)
    /// These conditions can be toggled on/off with CL flags
    pub filtered_by_flags: usize,

    /// The number of records that have been ignored because they failed the MAPQ filter.
    pub bad_mapq: usize,

    /// The number of records determined to be Paired-End.
    pub paired_end_reads: usize,

    /// The number of records determined to be Single-End.
    pub single_end_reads: usize,
}

/// General gene metrics that are tallied as a part of the
/// strandedness subcommand.
#[derive(Clone, Default, Serialize, Debug)]
pub struct GeneRecordMetrics {
    /// The total number of genes found in the GFF.
    pub total: usize,

    /// The number of genes that were found to be protein coding.
    /// If --all-genes is set this will not be tallied.
    pub protein_coding: usize,

    /// The number of genes which were discarded due to having
    /// an unknown or missing strand.
    pub bad_strand: usize,

    /// The number of genes randomly selected for evaluation.
    pub evaluated: usize,

    /// The number of genes which were discarded due to having
    /// mixed strands (the gene has exons on both strands).
    pub mixed_strands: usize,

    /// The number of genes which were discarded due to not having
    /// enough reads.
    pub not_enough_reads: usize,

    /// The number of genes considered for strandedness evidence.
    pub considered: usize,
}

/// General exon metrics that are tallied as a part of the
/// strandedness subcommand.
#[derive(Clone, Default, Serialize, Debug)]
pub struct ExonRecordMetrics {
    /// The total number of exons found in the GFF.
    pub total: usize,

    /// The number of exons discarded due to having an unknown/invalid strand.
    pub bad_strand: usize,
}

/// Struct for managing record tracking.
#[derive(Clone, Default, Debug)]
pub struct RecordTracker {
    /// Read metrics.
    pub reads: ReadRecordMetrics,

    /// Gene metrics.
    pub genes: GeneRecordMetrics,

    /// Exon metrics.
    pub exons: ExonRecordMetrics,
}

/// Struct holding the per read group results for an `ngs derive strandedness`
/// subcommand call.
#[derive(Debug, Serialize)]
pub struct ReadGroupDerivedStrandednessResult {
    /// Name of the read group.
    pub read_group: String,

    /// Whether or not strandedness was determined for this read group.
    pub succeeded: bool,

    /// The strandedness of this read group, if derivable.
    pub strandedness: Option<String>,

    /// The total number of reads in this read group.
    pub total: usize,

    /// The number of reads that are evidence of Forward Strandedness.
    pub forward: usize,

    /// The number of reads that are evidence of Reverse Strandedness.
    pub reverse: usize,

    /// The percent of evidence for Forward Strandedness.
    pub forward_pct: f64,

    /// The percent of evidence for Reverse Strandedness.
    pub reverse_pct: f64,
}

impl ReadGroupDerivedStrandednessResult {
    /// Creates a new [`ReadGroupDerivedStrandednessResult`].
    pub fn new(
        read_group: String,
        succeeded: bool,
        strandedness: Option<String>,
        forward: usize,
        reverse: usize,
    ) -> Self {
        ReadGroupDerivedStrandednessResult {
            read_group,
            succeeded,
            strandedness,
            total: forward + reverse,
            forward,
            reverse,
            forward_pct: (forward as f64 / (forward + reverse) as f64) * 100.0,
            reverse_pct: (reverse as f64 / (forward + reverse) as f64) * 100.0,
        }
    }
}

/// Struct holding the final results for an `ngs derive strandedness` subcommand
/// call.
#[derive(Debug, Serialize)]
pub struct DerivedStrandednessResult {
    /// Whether or not the `ngs derive strandedness` subcommand succeeded.
    pub succeeded: bool,

    /// The strandedness of this read group, if derivable.
    pub strandedness: Option<String>,

    /// The total number of reads.
    pub total: usize,

    /// The number of reads that are evidence of Forward Strandedness.
    pub forward: usize,

    /// The number of reads that are evidence of Reverse Strandedness.
    pub reverse: usize,

    /// The percent of evidence for Forward Strandedness.
    pub forward_pct: f64,

    /// The percent of evidence for Reverse Strandedness.
    pub reverse_pct: f64,

    /// Vector of [`ReadGroupDerivedStrandednessResult`]s.
    /// One for each read group in the BAM,
    /// and potentially one for any reads with an unknown read group.
    pub read_groups: Vec<ReadGroupDerivedStrandednessResult>,

    /// General read record metrics that are tallied as a part of the
    /// strandedness subcommand.
    pub read_metrics: ReadRecordMetrics,

    /// General gene metrics that are tallied as a part of the
    /// strandedness subcommand.
    pub gene_metrics: GeneRecordMetrics,

    /// General exon metrics that are tallied as a part of the
    /// strandedness subcommand.
    pub exon_metrics: ExonRecordMetrics,
}

impl DerivedStrandednessResult {
    /// Creates a new [`DerivedStrandednessResult`].
    pub fn new(
        succeeded: bool,
        strandedness: Option<String>,
        forward: usize,
        reverse: usize,
        read_groups: Vec<ReadGroupDerivedStrandednessResult>,
        metrics: RecordTracker,
    ) -> Self {
        DerivedStrandednessResult {
            succeeded,
            strandedness,
            total: forward + reverse,
            forward,
            reverse,
            forward_pct: (forward as f64 / (forward + reverse) as f64) * 100.0,
            reverse_pct: (reverse as f64 / (forward + reverse) as f64) * 100.0,
            read_groups,
            read_metrics: metrics.reads,
            gene_metrics: metrics.genes,
            exon_metrics: metrics.exons,
        }
    }
}
