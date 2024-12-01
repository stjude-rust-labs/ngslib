//! This module contains functions to validate the read group information in the header and the records.

use lazy_static::lazy_static;
use noodles::sam::alignment::Record;
use noodles::sam::header;
use noodles::sam::record::data::field::Tag;
use std::collections::HashSet;
use std::sync::Arc;
use tracing::warn;

/// Type alias for a read group pointer.
pub type ReadGroupPtr = Arc<String>;

lazy_static! {
    /// String used to represent an unknown read group. Wrapped in an Arc to prevent redundant memory usage.
    pub static ref UNKNOWN_READ_GROUP: ReadGroupPtr = Arc::new(String::from("unknown_read_group"));
}

/// Returns the read group tag from the record.
/// If the read group is not found in the record, the read group is set to "unknown_read_group".
/// TODO: Revisit this logic
pub fn get_read_group(
    record: &Record,
    found_rgs: Option<&mut HashSet<ReadGroupPtr>>,
) -> ReadGroupPtr {
    match (record.data().get(Tag::ReadGroup), found_rgs) {
        (Some(rg), Some(read_groups)) => {
            let rg = rg.to_string();
            if !read_groups.contains(&rg) {
                read_groups.insert(Arc::new(rg.clone()));
            }
            Arc::clone(read_groups.get(&rg).unwrap())
        }
        (Some(rg), None) => Arc::new(rg.to_string()),
        (None, _) => Arc::clone(&UNKNOWN_READ_GROUP),
    }
}

/// Compares the read group tags found in the records
/// and the read groups found in the header.
/// Returns a vector of read group names that were found in the header
/// but not in the records.
pub fn validate_read_group_info(
    found_rgs: &HashSet<ReadGroupPtr>,
    header: &header::Header,
) -> Vec<String> {
    let mut rgs_in_header_not_records = Vec::new();
    let mut rgs_in_records_not_header = Vec::new();

    for (rg_id, _) in header.read_groups() {
        if !found_rgs.contains(rg_id) {
            rgs_in_header_not_records.push(rg_id.to_string());
        }
    }
    if !rgs_in_header_not_records.is_empty() {
        warn!(
            "The following read groups were not found in the file: {:?}",
            rgs_in_header_not_records
        );
    }

    for rg_id in found_rgs {
        if !header.read_groups().contains_key(rg_id.as_str()) {
            rgs_in_records_not_header.push(rg_id.to_string());
        }
    }
    if !rgs_in_records_not_header.is_empty() {
        warn!(
            "The following read groups were not found in the header: {:?}",
            rgs_in_records_not_header
        );
    }

    rgs_in_header_not_records
}
