//! Functionality relating to the `ngs derive readlen` subcommand itself.

use anyhow::Context;
use clap::Args;
use num_format::{Locale, ToFormattedString};
use std::collections::HashMap;
use std::collections::HashSet;
use std::path::PathBuf;
use std::sync::Arc;
use tracing::info;

use crate::derive::readlen::compute;
use crate::utils::args::arg_in_range as cutoff_in_range;
use crate::utils::args::NumberOfRecords;
use crate::utils::display::RecordCounter;
use crate::utils::formats::bam::ParsedBAMFile;
use crate::utils::formats::utils::IndexCheck;
use crate::utils::read_groups::{get_read_group, validate_read_group_info, ReadGroupPtr};

/// Clap arguments for the `ngs derive readlen` subcommand.
#[derive(Args)]
pub struct DeriveReadlenArgs {
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

    /// Majority vote cutoff value as a fraction between [0.0, 1.0].
    #[arg(short, long, value_name = "F64", default_value = "0.7")]
    majority_vote_cutoff: f64,
}

/// Main function for the `ngs derive readlen` subcommand.
pub fn derive(args: DeriveReadlenArgs) -> anyhow::Result<()> {
    // (0) Parse arguments needed for subcommand.
    let majority_vote_cutoff = cutoff_in_range(args.majority_vote_cutoff, 0.0..=1.0)
        .with_context(|| "Majority vote cutoff is not within acceptable range")?;

    let mut read_lengths: HashMap<ReadGroupPtr, HashMap<usize, usize>> = HashMap::new();
    let mut found_rgs = HashSet::new();

    info!("Starting derive readlen subcommand.");

    let ParsedBAMFile {
        mut reader, header, ..
    } = crate::utils::formats::bam::open_and_parse(args.src, IndexCheck::None)?;

    // (1) Collect read lengths from reads within the
    // file. Support for sampling only a portion of the reads is provided.
    let mut counter = RecordCounter::default();
    for result in reader.records(&header.parsed) {
        let record = result?;
        let read_group = get_read_group(&record, Some(&mut found_rgs));
        let len = record.sequence().len();

        *read_lengths
            .entry(read_group)
            .or_default()
            .entry(len)
            .or_default() += 1;

        counter.inc();
        if counter.time_to_break(&args.num_records) {
            break;
        }
    }

    info!(
        "Processed {} records.",
        counter.get().to_formatted_string(&Locale::en)
    );

    // (2) Validate the read group information.
    let rgs_in_header_not_records = validate_read_group_info(&found_rgs, &header.parsed);
    for rg_id in rgs_in_header_not_records {
        read_lengths.insert(Arc::new(rg_id), HashMap::new());
    }

    // (3) Derive the consensus read length based on the read lengths gathered.
    let result = compute::predict(read_lengths, majority_vote_cutoff);

    // (4) Print the output to stdout as JSON (more support for different output
    // types may be added in the future, but for now, only JSON).
    let output = serde_json::to_string_pretty(&result).unwrap();
    println!("{}", output);

    anyhow::Ok(())
}
