//! Functionality relating to the `ngs derive instrument` subcommand itself.

use clap::Args;
use num_format::{Locale, ToFormattedString};
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::Arc;
use tracing::info;

use crate::derive::instrument::compute;
use crate::derive::instrument::reads::IlluminaReadName;
use crate::utils::args::NumberOfRecords;
use crate::utils::display::RecordCounter;
use crate::utils::formats::bam::ParsedBAMFile;
use crate::utils::formats::utils::IndexCheck;
use crate::utils::read_groups::{get_read_group, validate_read_group_info, ReadGroupPtr};

/// Clap arguments for the `ngs derive instrument` subcommand.
#[derive(Args)]
pub struct DeriveInstrumentArgs {
    // Source BAM.
    #[arg(value_name = "BAM")]
    src: PathBuf,

    /// Examine the first `n` records in the file.
    #[arg(
        short,
        long,
        default_value = "10000000",
        value_name = "'all' or a positive, non-zero integer"
    )]
    num_records: NumberOfRecords,
}

/// Main function for the `ngs derive instrument` subcommand.
pub fn derive(args: DeriveInstrumentArgs) -> anyhow::Result<()> {
    let src = args.src;
    let mut instrument_names: HashMap<ReadGroupPtr, HashSet<String>> = HashMap::new();
    let mut flowcell_names: HashMap<ReadGroupPtr, HashSet<String>> = HashMap::new();
    let mut metrics = compute::RecordMetrics::default();
    let mut found_rgs = HashSet::new();

    info!("Starting derive instrument subcommand.");

    let ParsedBAMFile {
        mut reader, header, ..
    } = crate::utils::formats::bam::open_and_parse(src, IndexCheck::None)?;

    // (1) Collect instrument names and flowcell names from reads within the
    // file. Support for sampling only a portion of the reads is provided.
    let mut counter = RecordCounter::default();
    for result in reader.records(&header.parsed) {
        let record = result?;
        let read_group = get_read_group(&record, Some(&mut found_rgs));

        if let Some(read_name) = record.read_name() {
            let name: &str = read_name.as_ref();

            match name.parse::<IlluminaReadName>() {
                Ok(read) => {
                    instrument_names
                        .entry(read_group.clone())
                        .or_default()
                        .insert(read.instrument_name);
                    metrics.found_instrument_name += 1;
                    // Get or init the flowcell set for this read group.
                    let fc_entry = flowcell_names.entry(read_group).or_default();
                    if let Some(fc) = read.flowcell {
                        fc_entry.insert(fc);
                        metrics.found_flowcell_name += 1;
                    }
                }
                Err(_) => {
                    metrics.bad_read_name += 1;
                }
            }
        } else {
            metrics.bad_read_name += 1;
        }

        counter.inc();
        if counter.time_to_break(&args.num_records) {
            break;
        }
    }

    info!(
        "Processed {} records.",
        counter.get().to_formatted_string(&Locale::en)
    );
    metrics.total_records = counter.get();

    // (2) Validate the read group information.
    let rgs_in_header_not_records = validate_read_group_info(&found_rgs, &header.parsed);
    for rg_id in rgs_in_header_not_records {
        let rg_ptr = Arc::new(rg_id);
        instrument_names.insert(rg_ptr.clone(), HashSet::new());
        flowcell_names.insert(rg_ptr, HashSet::new());
    }

    // (3) Derive the instrument results based on the detected
    // instrument names and flowcell names.
    let mut result = compute::predict(instrument_names, flowcell_names);
    result.records = metrics;

    // (4) Print the output to stdout as JSON (more support for different output
    // types may be added in the future, but for now, only JSON).
    let output = serde_json::to_string_pretty(&result).unwrap();
    println!("{}", output);

    Ok(())
}
