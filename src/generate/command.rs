//! Functionality related to the `ngs generate` subcommand itself.

use std::path::PathBuf;

use anyhow::Context;
use clap::ArgGroup;
use clap::Args;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use rand::prelude::*;
use tracing::info;

use crate::generate::providers::reference_provider::ReferenceGenomeSequenceProvider;
use crate::generate::providers::SequenceProvider;
use crate::utils::args::arg_in_range as error_rate_in_range;
use crate::utils::formats;

/// Command line arguments for `ngs generate`.
#[derive(Args)]
#[command(group(ArgGroup::new("record-count").required(true).args(["coverage", "num_records"])))]
pub struct GenerateArgs {
    /// Destination for the FASTQ file containing all read ones.
    read_ones_file: PathBuf,

    /// Destination for the FASTQ file containing all read twos.
    read_twos_file: PathBuf,

    /// One or more reference FASTAs to generate the data based off of.
    #[arg(required = true)] // required implies one or more
    reference_providers: Vec<String>,

    /// The error rate for the sequencer as a fraction between [0.0, 1.0] (per base).
    #[arg(short, long, value_name = "F64", default_value = "0.0001")]
    error_rate: f64,

    /// Specifies the number of records to generate.
    #[arg(short, long, value_name = "USIZE", conflicts_with = "coverage")]
    num_records: Option<usize>,

    /// Dynamically calculate the number of reads needed for a particular mean coverage.
    #[arg(short, long, value_name = "USIZE", conflicts_with = "num_records")]
    coverage: Option<usize>,
}

/// Main function for the `ngs generate` subcommand.
pub fn generate(args: GenerateArgs) -> anyhow::Result<()> {
    // (0) Parse arguments needed for subcommand.
    let _error_rate = error_rate_in_range(args.error_rate, 0.0..=1.0)
        .with_context(|| "Error rate is not within acceptable range")?;

    let result: anyhow::Result<Vec<_>> = args
        .reference_providers
        .iter()
        .map(|provider_as_string| provider_as_string.parse::<ReferenceGenomeSequenceProvider>())
        .collect();
    let reference_providers = result.with_context(|| "parsing reference providers")?;

    let reads_one_file = args.read_ones_file;
    let reads_two_file = args.read_twos_file;

    info!("Starting generate command...");
    let mut writer_read_one = formats::fastq::writer(&reads_one_file)
        .with_context(|| format!("opening reads one file: {}", reads_one_file.display()))?;
    let mut writer_read_two = formats::fastq::writer(&reads_two_file)
        .with_context(|| format!("opening reads two file: {}", reads_two_file.display()))?;

    let mut total_reads: usize = 0;

    if let Some(num_records) = args.num_records {
        total_reads = num_records
    } else if let Some(coverage) = args.coverage {
        total_reads = reference_providers[0].reads_needed_for_coverage(coverage);
    }

    // (2) Set up the output writers.
    info!("Generating {} reads...", total_reads);

    // (3) Set up the progress bar.
    let pb = ProgressBar::new(total_reads as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{prefix:.cyan.bold} {spinner:.green} [{elapsed_precise}] [{bar}] {pos}/{len} ({per_sec}, {eta})")
            .progress_chars("=> "),
    );
    pb.set_prefix("Generating");

    // (4) Generate the reads and write them to their respective files.
    let mut rng = ThreadRng::default();

    let mut i = 0;
    while i < total_reads {
        let selected_genome = reference_providers
            .choose_weighted(&mut rng, |x| x.weight)
            .unwrap();
        let read_pair = selected_genome.generate_read_pair(
            format!("ngs:{}", selected_genome.filename),
            (i + 1).to_string(),
        );
        writer_read_one
            .write_record(read_pair.get_forward_read())
            .with_context(|| "could not write record to read one file")?;
        writer_read_two
            .write_record(read_pair.get_reverse_read())
            .with_context(|| "could not write record to read two file")?;

        if i > 0 && i % 5_000 == 0 {
            pb.inc(5000);
        }

        i += 1;
    }

    pb.set_style(
        ProgressStyle::default_bar()
            .template("{prefix:.green.bold} {msg:.white.bold} [{elapsed_precise}] [{bar}] {pos}/{len} ({per_sec}, {eta})")
            .progress_chars("=> "),
    );
    pb.set_prefix("✓");
    pb.finish_with_message("Finished");
    Ok(())
}
