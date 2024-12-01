//! Module holding the results structs for the `ngs derive endedness` subcommand.

use serde::Serialize;

use crate::derive::endedness::compute::OrderingFlagsCounts;

/// Struct holding the per read group results for an `ngs derive endedness`
/// subcommand call.
#[derive(Debug, Serialize)]
pub struct ReadGroupDerivedEndednessResult {
    /// Name of the read group.
    pub read_group: String,

    /// Whether or not an endedness was determined for this read group.
    pub succeeded: bool,

    /// The endedness of this read group, if derivable.
    pub endedness: Option<String>,

    /// The number of reads without 0x1 set.
    pub unsegmented: usize,

    /// The f+l- read count.
    pub first: usize,

    /// The f-l+ read count.
    pub last: usize,

    /// The f+l+ read count.
    pub both: usize,

    /// The f-l- read count.
    pub neither: usize,

    /// The reads per template (RPT).
    /// Only available if `args.calc_rpt` is true.
    pub rpt: Option<f64>,
}

impl ReadGroupDerivedEndednessResult {
    /// Creates a new [`ReadGroupDerivedEndednessResult`].
    pub fn new(
        read_group: String,
        succeeded: bool,
        endedness: Option<String>,
        counts: OrderingFlagsCounts,
        rpt: Option<f64>,
    ) -> Self {
        ReadGroupDerivedEndednessResult {
            read_group,
            succeeded,
            endedness,
            unsegmented: counts.unsegmented,
            first: counts.first,
            last: counts.last,
            both: counts.both,
            neither: counts.neither,
            rpt,
        }
    }
}

/// Struct holding the final results for an `ngs derive endedness` subcommand
/// call.
#[derive(Debug, Serialize)]
pub struct DerivedEndednessResult {
    /// Whether or not the `ngs derive endedness` subcommand succeeded.
    pub succeeded: bool,

    /// The overall endedness, if derivable.
    pub endedness: Option<String>,

    /// The number of reads without 0x1 set.
    pub unsegmented: usize,

    /// The overall f+l- read count.
    pub first: usize,

    /// The overall f-l+ read count.
    pub last: usize,

    /// The overall f+l+ read count.
    pub both: usize,

    /// The overall f-l- read count.
    pub neither: usize,

    /// The overall reads per template (RPT).
    /// Only available if `args.calc_rpt` is true.
    pub rpt: Option<f64>,

    /// Vector of [`ReadGroupDerivedEndednessResult`]s.
    /// One for each read group in the BAM,
    /// and potentially one for any reads with an unknown read group.
    pub read_groups: Vec<ReadGroupDerivedEndednessResult>,
}

impl DerivedEndednessResult {
    /// Creates a new [`DerivedEndednessResult`].
    pub fn new(
        succeeded: bool,
        endedness: Option<String>,
        counts: OrderingFlagsCounts,
        rpt: Option<f64>,
        read_groups: Vec<ReadGroupDerivedEndednessResult>,
    ) -> Self {
        DerivedEndednessResult {
            succeeded,
            endedness,
            unsegmented: counts.unsegmented,
            first: counts.first,
            last: counts.last,
            both: counts.both,
            neither: counts.neither,
            rpt,
            read_groups,
        }
    }
}
