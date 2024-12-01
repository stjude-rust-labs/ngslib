//! Module holding the logic for computing the endedness of a BAM.

use std::collections::{HashMap, HashSet};
use std::ops::{Add, AddAssign};
use std::sync::Arc;
use tracing::warn;

use crate::derive::endedness::results;
use crate::utils::read_groups::ReadGroupPtr;

/// Struct holding the ordering flags for a single read group.
#[derive(Debug, Clone, Default)]
pub struct OrderingFlagsCounts {
    /// The number of reads without 0x1 set.
    pub unsegmented: usize,

    /// The number of reads with the first in template flag set.
    pub first: usize,

    /// The number of reads with the last in template flag set.
    pub last: usize,

    /// The number of reads with both the first and last in template flags set.
    pub both: usize,

    /// The number of reads with neither the first nor last in template flags set.
    pub neither: usize,
}
impl OrderingFlagsCounts {
    /// Creates a new [`OrderingFlagsCounts`].
    pub fn new() -> Self {
        OrderingFlagsCounts {
            unsegmented: 0,
            first: 0,
            last: 0,
            both: 0,
            neither: 0,
        }
    }
}

impl Add for OrderingFlagsCounts {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        OrderingFlagsCounts {
            unsegmented: self.unsegmented + other.unsegmented,
            first: self.first + other.first,
            last: self.last + other.last,
            both: self.both + other.both,
            neither: self.neither + other.neither,
        }
    }
}

impl AddAssign for OrderingFlagsCounts {
    fn add_assign(&mut self, other: Self) {
        self.unsegmented += other.unsegmented;
        self.first += other.first;
        self.last += other.last;
        self.both += other.both;
        self.neither += other.neither;
    }
}

/// Calculate the reads per template overall and for each read group.
fn calculate_reads_per_template(
    read_names: HashMap<String, Vec<ReadGroupPtr>>,
    reads_per_template: &mut HashMap<ReadGroupPtr, f64>,
) -> f64 {
    let mut total_reads: usize = 0;
    let mut total_templates: usize = 0;
    let mut read_group_reads: HashMap<ReadGroupPtr, usize> = HashMap::new();
    let mut read_group_templates: HashMap<ReadGroupPtr, usize> = HashMap::new();

    let mut warning_count: usize = 0;

    for (read_name, read_groups) in read_names.iter() {
        let num_reads = read_groups.len();
        total_reads += num_reads;
        total_templates += 1;

        let read_group_set: HashSet<ReadGroupPtr> = read_groups.iter().cloned().collect();

        if read_group_set.len() == 1 {
            // All found read groups assigned to this QNAME are the same.
            // We assume this means all the reads came from the same template.
            let read_group = Arc::clone(read_group_set.iter().next().unwrap());

            *read_group_reads.entry(Arc::clone(&read_group)).or_default() += num_reads;
            *read_group_templates.entry(read_group).or_default() += 1;
        } else {
            // The QNAME is in multiple read groups.
            // We assume this means the reads came from multiple templates.
            // More specifically, we assume that exactly one template will originate from each read group.
            warning_count += 1;
            match warning_count {
                1..=100 => {
                    warn!(
                        "QNAME: '{}' is in multiple read groups: {:?}",
                        read_name, read_group_set
                    );
                }
                101 => warn!(
                    "Too many warnings about QNAMEs in multiple read groups. Stopping warnings."
                ),
                _ => (),
            }

            for read_group in read_groups {
                *read_group_reads.entry(Arc::clone(read_group)).or_default() += 1;
            }
            for read_group in read_group_set {
                *read_group_templates.entry(read_group).or_default() += 1;
            }
        }
    }

    if warning_count > 100 {
        warn!(
            "{} QNAMEs were found in multiple read groups.",
            warning_count
        );
    }

    let overall_rpt = total_reads as f64 / total_templates as f64;

    for (read_group, num_reads) in read_group_reads.iter() {
        let num_templates = read_group_templates.get(read_group).unwrap();
        let rpt = *num_reads as f64 / *num_templates as f64;
        reads_per_template.insert(Arc::clone(read_group), rpt);
    }

    overall_rpt
}

fn predict_endedness(
    read_group_name: String,
    rg_ordering_flags: &OrderingFlagsCounts,
    paired_deviance: f64,
    reads_per_template: Option<f64>,
    round_rpt: bool,
) -> results::ReadGroupDerivedEndednessResult {
    let unsegmented = rg_ordering_flags.unsegmented;
    let first = rg_ordering_flags.first;
    let last = rg_ordering_flags.last;
    let both = rg_ordering_flags.both;
    let neither = rg_ordering_flags.neither;

    // all zeroes (Perform this check before creating the result struct
    // so that we don't have to clone the read group name)
    if unsegmented == 0 && first == 0 && last == 0 && both == 0 && neither == 0 {
        warn!(
            "No reads were detected in this read group: {}",
            read_group_name
        );
        return results::ReadGroupDerivedEndednessResult::new(
            read_group_name,
            false,
            None,
            rg_ordering_flags.clone(),
            reads_per_template,
        );
    }

    let mut result = results::ReadGroupDerivedEndednessResult::new(
        read_group_name,
        false,
        None,
        rg_ordering_flags.clone(),
        reads_per_template,
    );

    // only unsegmented present
    if unsegmented > 0 && first == 0 && last == 0 && both == 0 && neither == 0 {
        match reads_per_template {
            Some(rpt) => {
                if rpt == 1.0 || (round_rpt && rpt.round() as usize == 1) {
                    result.succeeded = true;
                    result.endedness = Some(String::from("Single-End"));
                }
            }
            None => {
                result.succeeded = true;
                result.endedness = Some(String::from("Single-End"));
            }
        }
        return result;
    }
    // unsegmented reads are present, and so are other types of reads.
    if unsegmented > 0 {
        return result;
    }
    // now unsegmented is guarenteed to be 0

    // only first present
    if first > 0 && last == 0 && both == 0 && neither == 0 {
        return result;
    }
    // only last present
    if first == 0 && last > 0 && both == 0 && neither == 0 {
        return result;
    }
    // only both present
    if first == 0 && last == 0 && both > 0 && neither == 0 {
        // Prior logic (before addition of unsegmented checks) left as comment for posterity
        // match reads_per_template {
        //     Some(rpt) => {
        //         if rpt == 1.0 || (round_rpt && rpt.round() as usize == 1) {
        //             result.succeeded = true;
        //             result.endedness = String::from("Single-End");
        //         }
        //     }
        //     None => {
        //         result.succeeded = true;
        //         result.endedness = String::from("Single-End");
        //     }
        // }
        return result;
    }
    // only neither present
    if first == 0 && last == 0 && both == 0 && neither > 0 {
        return result;
    }
    // first/last mixed with both/neither
    if (first > 0 || last > 0) && (both > 0 || neither > 0) {
        return result;
    }
    // any mix of both/neither, regardless of first/last
    if both > 0 && neither > 0 {
        return result;
    }

    // both and neither are now guarenteed to be 0
    // We only need to check first and last

    let first_frac = first as f64 / (first + last) as f64;
    let lower_limit = 0.5 - paired_deviance;
    let upper_limit = 0.5 + paired_deviance;
    if (first == last) || (lower_limit <= first_frac && first_frac <= upper_limit) {
        match reads_per_template {
            Some(rpt) => {
                if rpt == 2.0 || (round_rpt && rpt.round() as usize == 2) {
                    result.succeeded = true;
                    result.endedness = Some(String::from("Paired-End"));
                }
            }
            None => {
                result.succeeded = true;
                result.endedness = Some(String::from("Paired-End"));
            }
        }
    }
    result
}

/// Main method to evaluate the collected ordering flags and
/// return a result for the endedness of the file. This may fail, and the
/// resulting [`results::DerivedEndednessResult`] should be evaluated accordingly.
pub fn predict(
    ordering_flags: HashMap<ReadGroupPtr, OrderingFlagsCounts>,
    read_names: Option<HashMap<String, Vec<ReadGroupPtr>>>,
    paired_deviance: f64,
    round_rpt: bool,
) -> results::DerivedEndednessResult {
    let mut rg_rpts: HashMap<ReadGroupPtr, f64> = HashMap::new();
    let mut overall_rpt: Option<f64> = None;
    if let Some(read_names) = read_names {
        overall_rpt = Some(calculate_reads_per_template(read_names, &mut rg_rpts));
    }

    let mut overall_flags = OrderingFlagsCounts::new();
    let mut rg_results = Vec::new();

    for (read_group, rg_ordering_flags) in ordering_flags.iter() {
        overall_flags += rg_ordering_flags.clone();

        let result = predict_endedness(
            read_group.to_string(),
            rg_ordering_flags,
            paired_deviance,
            rg_rpts.get(read_group).copied(),
            round_rpt,
        );
        rg_results.push(result);
    }

    let overall_result = predict_endedness(
        "overall".to_string(),
        &overall_flags,
        paired_deviance,
        overall_rpt,
        round_rpt,
    );

    results::DerivedEndednessResult::new(
        overall_result.succeeded,
        overall_result.endedness,
        overall_flags,
        overall_rpt,
        rg_results,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_predict_endedness_from_first_and_last() {
        let mut ordering_flags: HashMap<ReadGroupPtr, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(
            Arc::new("overall".to_string()),
            OrderingFlagsCounts {
                unsegmented: 0,
                first: 1,
                last: 1,
                both: 0,
                neither: 0,
            },
        );
        let result = predict_endedness(
            "overall".to_string(),
            ordering_flags
                .get(&Arc::new("overall".to_string()))
                .unwrap(),
            0.0,
            None,
            false,
        );
        assert!(result.succeeded);
        assert_eq!(result.endedness, Some("Paired-End".to_string()));
        assert_eq!(result.first, 1);
        assert_eq!(result.last, 1);
        assert_eq!(result.both, 0);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, None);
    }

    #[test]
    fn test_predict_endedness_from_unsegmented() {
        let mut ordering_flags: HashMap<ReadGroupPtr, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(
            Arc::new("overall".to_string()),
            OrderingFlagsCounts {
                unsegmented: 1,
                first: 0,
                last: 0,
                both: 0,
                neither: 0,
            },
        );
        let result = predict_endedness(
            "overall".to_string(),
            ordering_flags
                .get(&Arc::new("overall".to_string()))
                .unwrap(),
            0.0,
            None,
            false,
        );
        assert!(result.succeeded);
        assert_eq!(result.endedness, Some("Single-End".to_string()));
        assert_eq!(result.first, 0);
        assert_eq!(result.last, 0);
        assert_eq!(result.both, 0);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, None);
    }

    #[test]
    fn test_predict_endedness_from_all_zero_counts() {
        let mut ordering_flags: HashMap<ReadGroupPtr, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(Arc::new(String::from("rg1")), OrderingFlagsCounts::new());
        let result = predict_endedness(
            String::from("rg1"),
            ordering_flags.get(&Arc::new(String::from("rg1"))).unwrap(),
            0.0,
            None,
            false,
        );
        assert!(!result.succeeded);
        assert_eq!(result.endedness, None);
        assert_eq!(result.first, 0);
        assert_eq!(result.last, 0);
        assert_eq!(result.both, 0);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, None);
    }

    #[test]
    fn test_predict_from_only_first() {
        let mut ordering_flags: HashMap<ReadGroupPtr, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(
            Arc::new("overall".to_string()),
            OrderingFlagsCounts {
                unsegmented: 0,
                first: 1,
                last: 0,
                both: 0,
                neither: 0,
            },
        );
        let result = predict(ordering_flags, None, 0.0, false);
        assert!(!result.succeeded);
        assert_eq!(result.endedness, None);
        assert_eq!(result.first, 1);
        assert_eq!(result.last, 0);
        assert_eq!(result.both, 0);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, None);
        assert_eq!(result.read_groups.len(), 1);
    }

    #[test]
    fn test_predict_from_only_last() {
        let mut ordering_flags: HashMap<ReadGroupPtr, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(
            Arc::new("overall".to_string()),
            OrderingFlagsCounts {
                unsegmented: 0,
                first: 0,
                last: 1,
                both: 0,
                neither: 0,
            },
        );
        let result = predict(ordering_flags, None, 0.0, false);
        assert!(!result.succeeded);
        assert_eq!(result.endedness, None);
        assert_eq!(result.first, 0);
        assert_eq!(result.last, 1);
        assert_eq!(result.both, 0);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, None);
        assert_eq!(result.read_groups.len(), 1);
    }

    #[test]
    fn test_predict_from_only_both() {
        let mut ordering_flags: HashMap<ReadGroupPtr, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(
            Arc::new("overall".to_string()),
            OrderingFlagsCounts {
                unsegmented: 0,
                first: 0,
                last: 0,
                both: 1,
                neither: 0,
            },
        );
        let result = predict(ordering_flags, None, 0.0, false);
        assert!(!result.succeeded);
        assert_eq!(result.endedness, None);
        assert_eq!(result.first, 0);
        assert_eq!(result.last, 0);
        assert_eq!(result.both, 1);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, None);
        assert_eq!(result.read_groups.len(), 1);
    }

    #[test]
    fn test_predict_from_only_neither() {
        let mut ordering_flags: HashMap<ReadGroupPtr, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(
            Arc::new("overall".to_string()),
            OrderingFlagsCounts {
                unsegmented: 0,
                first: 0,
                last: 0,
                both: 0,
                neither: 1,
            },
        );
        let result = predict(ordering_flags, None, 0.0, false);
        assert!(!result.succeeded);
        assert_eq!(result.endedness, None);
        assert_eq!(result.first, 0);
        assert_eq!(result.last, 0);
        assert_eq!(result.both, 0);
        assert_eq!(result.neither, 1);
        assert_eq!(result.rpt, None);
        assert_eq!(result.read_groups.len(), 1);
    }

    #[test]
    fn test_predict_from_first_and_last() {
        let mut ordering_flags: HashMap<ReadGroupPtr, OrderingFlagsCounts> = HashMap::new();
        ordering_flags.insert(
            Arc::new("overall".to_string()),
            OrderingFlagsCounts {
                unsegmented: 0,
                first: 1,
                last: 1,
                both: 0,
                neither: 0,
            },
        );
        let result = predict(ordering_flags, None, 0.0, false);
        assert!(result.succeeded);
        assert_eq!(result.endedness, Some("Paired-End".to_string()));
        assert_eq!(result.first, 1);
        assert_eq!(result.last, 1);
        assert_eq!(result.both, 0);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, None);
        assert_eq!(result.read_groups.len(), 1);
    }

    #[test]
    fn test_calculate_reads_per_template() {
        let mut read_names: HashMap<String, Vec<ReadGroupPtr>> = HashMap::new();
        let rg_paired = Arc::new("rg_paired".to_string());
        let rg_single = Arc::new("rg_single".to_string());
        read_names.insert(
            "read1".to_string(),
            vec![Arc::clone(&rg_paired), Arc::clone(&rg_paired)],
        );
        read_names.insert(
            "read2".to_string(),
            vec![
                Arc::clone(&rg_paired),
                Arc::clone(&rg_paired),
                Arc::clone(&rg_single),
            ],
        );
        read_names.insert("read3".to_string(), vec![Arc::clone(&rg_single)]);
        read_names.insert(
            "read4".to_string(),
            vec![Arc::clone(&rg_paired), Arc::clone(&rg_paired)],
        );
        read_names.insert(
            "read5".to_string(),
            vec![
                Arc::clone(&rg_paired),
                Arc::clone(&rg_paired),
                Arc::clone(&rg_single),
            ],
        );
        let mut rg_rpts: HashMap<ReadGroupPtr, f64> = HashMap::new();
        let overall_rpt = calculate_reads_per_template(read_names, &mut rg_rpts);
        assert_eq!(rg_rpts.len(), 2);
        assert_eq!(overall_rpt, 2.2);
        assert_eq!(rg_rpts.get(&Arc::clone(&rg_paired)).unwrap(), &2.0);
        assert_eq!(rg_rpts.get(&Arc::clone(&rg_single)).unwrap(), &1.0);
    }

    #[test]
    fn test_predict_with_rpt_complex() {
        let mut ordering_flags: HashMap<ReadGroupPtr, OrderingFlagsCounts> = HashMap::new();
        let rg_paired = Arc::new("rg_paired".to_string());
        let rg_single = Arc::new("rg_single".to_string());
        ordering_flags.insert(
            Arc::clone(&rg_paired),
            OrderingFlagsCounts {
                unsegmented: 0,
                first: 8,
                last: 8,
                both: 0,
                neither: 0,
            },
        );
        ordering_flags.insert(
            Arc::clone(&rg_single),
            OrderingFlagsCounts {
                unsegmented: 2,
                first: 0,
                last: 0,
                both: 0,
                neither: 0,
            },
        );
        let mut read_names: HashMap<String, Vec<ReadGroupPtr>> = HashMap::new();
        read_names.insert(
            "read1".to_string(),
            vec![Arc::clone(&rg_paired), Arc::clone(&rg_paired)],
        );
        read_names.insert(
            "read2".to_string(),
            vec![
                Arc::clone(&rg_paired),
                Arc::clone(&rg_paired),
                Arc::clone(&rg_single),
            ],
        );
        read_names.insert("read3".to_string(), vec![Arc::clone(&rg_single)]);
        read_names.insert(
            "read4".to_string(),
            vec![Arc::clone(&rg_paired), Arc::clone(&rg_paired)],
        );
        read_names.insert(
            "read5".to_string(),
            vec![
                Arc::clone(&rg_paired),
                Arc::clone(&rg_paired),
                Arc::clone(&rg_single),
            ],
        );
        let result = predict(ordering_flags, Some(read_names), 0.0, false);
        assert!(!result.succeeded);
        assert_eq!(result.endedness, None);
        assert_eq!(result.unsegmented, 2);
        assert_eq!(result.first, 8);
        assert_eq!(result.last, 8);
        assert_eq!(result.both, 0);
        assert_eq!(result.neither, 0);
        assert_eq!(result.rpt, Some(2.2));
        assert_eq!(result.read_groups.len(), 2);
        // We can't know which read group will be first in the vector.
        // But both should succeed.
        assert!(result.read_groups[0].succeeded && result.read_groups[1].succeeded);
    }
}
