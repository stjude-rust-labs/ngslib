//! Module holding the logic for computing the consensus read length.

use serde::Serialize;
use std::collections::HashMap;
use tracing::warn;

use crate::utils::read_groups::ReadGroupPtr;

/// Struct holding the per read group results for an `ngs derive readlen`
/// subcommand call.
#[derive(Debug, Serialize)]
pub struct ReadGroupDerivedReadlenResult {
    /// The read group that these results are associated with.
    pub read_group: String,

    /// Whether or not the `ngs derive readlen` subcommand succeeded
    /// for this read group.
    pub succeeded: bool,

    /// The consensus read length, if derivable.
    pub consensus_read_length: Option<usize>,

    /// The majority vote percentage of the consensus read length.
    pub majority_pct_detected: f64,

    /// Status of the evidence that supports (or does not support) the
    /// consensus read length.
    pub evidence: Vec<(usize, usize)>,
}

impl ReadGroupDerivedReadlenResult {
    /// Creates a new [`ReadGroupDerivedReadlenResult`].
    pub fn new(
        read_group: String,
        succeeded: bool,
        consensus_read_length: Option<usize>,
        majority_pct_detected: f64,
        evidence: Vec<(usize, usize)>,
    ) -> Self {
        ReadGroupDerivedReadlenResult {
            read_group,
            succeeded,
            consensus_read_length,
            majority_pct_detected,
            evidence,
        }
    }
}

/// Struct holding the final results for an `ngs derive readlen` subcommand
/// call.
#[derive(Debug, Serialize)]
pub struct DerivedReadlenResult {
    /// Whether or not the `ngs derive readlen` subcommand succeeded.
    pub succeeded: bool,

    /// The consensus read length, if derivable.
    pub consensus_read_length: Option<usize>,

    /// The majority vote percentage of the consensus read length.
    pub majority_pct_detected: f64,

    /// Vector of [`ReadGroupDerivedReadlenResult`]s.
    /// One for each read group in the BAM,
    /// and potentially one for any reads with an unknown read group.
    pub read_groups: Vec<ReadGroupDerivedReadlenResult>,

    /// Status of the evidence that supports (or does not support) the
    /// consensus read length.
    pub evidence: Vec<(usize, usize)>,
}

impl DerivedReadlenResult {
    /// Creates a new [`DerivedReadlenResult`].
    pub fn new(
        succeeded: bool,
        consensus_read_length: Option<usize>,
        majority_pct_detected: f64,
        read_groups: Vec<ReadGroupDerivedReadlenResult>,
        evidence: Vec<(usize, usize)>,
    ) -> Self {
        DerivedReadlenResult {
            succeeded,
            consensus_read_length,
            majority_pct_detected,
            read_groups,
            evidence,
        }
    }
}

/// Predicts the consensus read length for a given read group based on the
/// read lengths and a majority vote cutoff.
pub fn predict_readlen(
    read_group: String,
    read_lengths: &HashMap<usize, usize>,
    majority_vote_cutoff: f64,
) -> ReadGroupDerivedReadlenResult {
    let mut read_lengths: Vec<(usize, usize)> =
        read_lengths.iter().map(|(k, v)| (*k, *v)).collect();

    read_lengths.sort_by(|a, b| b.0.cmp(&a.0));

    // Tally the number of reads
    let num_reads: usize = read_lengths.iter().map(|(_, count)| count).sum();

    let (majority_detected, consensus_read_length) = match num_reads == 0 {
        true => {
            warn!("No reads were detected for read group: {}", read_group);
            (0.0, None)
        }
        false => (
            read_lengths[0].1 as f64 / num_reads as f64,
            Some(read_lengths[0].0),
        ),
    };

    match majority_detected >= majority_vote_cutoff {
        true => ReadGroupDerivedReadlenResult::new(
            read_group,
            true,
            consensus_read_length,
            majority_detected * 100.0,
            read_lengths,
        ),
        false => ReadGroupDerivedReadlenResult::new(
            read_group,
            false,
            None,
            majority_detected * 100.0,
            read_lengths,
        ),
    }
}

/// Main method to evaluate the collected read lengths and
/// return a result for the consensus read length. This may fail, and the
/// resulting [`DerivedReadlenResult`] should be evaluated accordingly.
pub fn predict(
    read_lengths: HashMap<ReadGroupPtr, HashMap<usize, usize>>,
    majority_vote_cutoff: f64,
) -> DerivedReadlenResult {
    // Iterate over the read lengths and predict the consensus read length.
    let mut rg_results = Vec::new();
    let mut overall_lengths = HashMap::new();

    for (read_group, lengths) in read_lengths {
        let result = predict_readlen(read_group.to_string(), &lengths, majority_vote_cutoff);
        rg_results.push(result);

        for (length, count) in lengths {
            *overall_lengths.entry(length).or_default() += count;
        }
    }

    let overall_result = predict_readlen(
        "overall".to_string(),
        &overall_lengths,
        majority_vote_cutoff,
    );

    // Sort the read lengths by their key for output.
    let mut overall_lengths: Vec<(usize, usize)> = overall_lengths.into_iter().collect();
    overall_lengths.sort_by(|a, b| b.0.cmp(&a.0));

    DerivedReadlenResult::new(
        overall_result.succeeded,
        overall_result.consensus_read_length,
        overall_result.majority_pct_detected,
        rg_results,
        overall_lengths,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn test_derive_readlen_from_empty_hashmap() {
        let read_lengths = HashMap::new();
        let result = predict(read_lengths, 0.7);
        assert!(!result.succeeded);
        assert_eq!(result.consensus_read_length, None);
        assert_eq!(result.majority_pct_detected, 0.0);
        assert_eq!(result.evidence, Vec::new());
    }

    #[test]
    fn test_derive_readlen_when_all_readlengths_equal() {
        let read_lengths =
            HashMap::from([(Arc::new("RG1".to_string()), HashMap::from([(100, 10)]))]);
        let result = predict(read_lengths, 1.0);
        assert!(result.succeeded);
        assert_eq!(result.consensus_read_length, Some(100));
        assert_eq!(result.majority_pct_detected, 100.0);
        assert_eq!(result.evidence, Vec::from([(100, 10)]));
    }

    #[test]
    fn test_derive_readlen_success_when_not_all_readlengths_equal() {
        let read_lengths = HashMap::from([(
            Arc::new("RG1".to_string()),
            HashMap::from([(101, 1000), (100, 5), (99, 5)]),
        )]);
        let result = predict(read_lengths, 0.7);
        assert!(result.succeeded);
        assert_eq!(result.consensus_read_length, Some(101));
        assert!(result.majority_pct_detected > 99.0);
        assert_eq!(result.evidence, Vec::from([(101, 1000), (100, 5), (99, 5)]));
    }

    #[test]
    fn test_derive_readlen_fail_when_not_all_readlengths_equal() {
        let read_lengths = HashMap::from([
            (
                Arc::new("RG1".to_string()),
                HashMap::from([(101, 5), (99, 5)]),
            ),
            (Arc::new("RG2".to_string()), HashMap::from([(100, 1000)])),
        ]);
        let result = predict(read_lengths, 0.7);
        assert!(!result.succeeded);
        assert_eq!(result.consensus_read_length, None);
        assert!(result.majority_pct_detected < 0.7);
        assert_eq!(result.evidence, Vec::from([(101, 5), (100, 1000), (99, 5)]));
    }
}
