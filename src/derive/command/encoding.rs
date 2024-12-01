//! Functionality relating to the `ngs derive encoding` subcommand itself.

use anyhow::{Context, Ok};
use clap::Args;
use noodles::bam;
use num_format::{Locale, ToFormattedString};
use std::collections::HashSet;
use std::io::BufReader;
use std::path::PathBuf;
use tracing::info;

use crate::derive::encoding::compute;
use crate::utils::args::NumberOfRecords;
use crate::utils::display::RecordCounter;

/// Clap arguments for the `ngs derive encoding` subcommand.
#[derive(Args)]
pub struct DeriveEncodingArgs {
    /// Source BAM.
    #[arg(value_name = "BAM")]
    src: PathBuf,

    /// Examine the first `n` records in the file.
    #[arg(
        short,
        long,
        default_value_t,
        value_name = "'all' or a positive, non-zero integer"
    )]
    num_records: NumberOfRecords,
}

/// Main function for the `ngs derive encoding` subcommand.
pub fn derive(args: DeriveEncodingArgs) -> anyhow::Result<()> {
    info!("Starting derive encoding subcommand.");

    let file = std::fs::File::open(args.src);
    let reader = file
        .map(BufReader::new)
        .with_context(|| "opening BAM file")?;
    let mut reader = bam::Reader::new(reader);
    let _header: String = reader.read_header()?.parse()?;
    reader.read_reference_sequences()?;

    let mut score_set: HashSet<u8> = HashSet::new();

    // (1) Collect quality scores from reads within the
    // file. Support for sampling only a portion of the reads is provided.
    let mut counter = RecordCounter::default();
    for result in reader.lazy_records() {
        let record = result?;

        for i in 0..record.quality_scores().len() {
            let score = record.quality_scores().as_ref()[i];
            score_set.insert(score);
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

    // (2) Derive encoding from the observed quality scores
    let result = compute::predict(score_set)?;

    // (3) Print the output to stdout as JSON (more support for different output
    // types may be added in the future, but for now, only JSON).
    let output = serde_json::to_string_pretty(&result).unwrap();
    println!("{}", output);

    Ok(())
}
