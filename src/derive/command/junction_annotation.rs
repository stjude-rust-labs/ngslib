//! Functionality relating to the `ngs derive junction-annotation` subcommand itself.

use anyhow::Context;
use clap::Args;
use noodles::sam::record::MappingQuality;
use num_format::{Locale, ToFormattedString};
use std::collections::HashMap;
use std::path::PathBuf;
use tracing::{debug, info};

use crate::derive::junction_annotation::compute;
use crate::derive::junction_annotation::results::JunctionAnnotationResults;
use crate::utils::display::RecordCounter;
use crate::utils::formats;
use crate::utils::formats::bam::ParsedBAMFile;
use crate::utils::formats::utils::IndexCheck;

/// Clap arguments for the `ngs derive junction-annotation` subcommand.
#[derive(Args)]
pub struct JunctionAnnotationArgs {
    /// Source BAM.
    #[arg(value_name = "BAM")]
    src: PathBuf,

    /// Features GFF file.
    #[arg(short = 'f', long, required = true, value_name = "PATH")]
    features_gff: PathBuf,

    /// Name of the exon region feature for the gene model used.
    #[arg(long, value_name = "STRING", default_value = "exon")]
    exon_feature_name: String,

    /// Minimum intron length to consider.
    /// An intron is defined as an `N` CIGAR operation of any length.
    #[arg(short = 'i', long, value_name = "USIZE", default_value = "50")]
    min_intron_length: usize,

    /// Minimum number of reads supporting a junction to be considered.
    #[arg(short = 'r', long, value_name = "USIZE", default_value = "2")]
    min_read_support: usize,

    /// Minumum mapping quality for a record to be considered.
    /// Default behavior is to ignore MAPQ values,
    /// which allows reads with _missing_ MAPQs to be considered.
    /// Specify any u8 value (lower than 255) to enable this filter.
    /// Some aligners erroneously use 255 as the score for a uniquely mapped read;
    /// however, 255 is reserved by the spec for a missing MAPQ value.
    /// Therefore BAMs produced by aligners using 255 erroneously
    /// are not compatible with setting this option.
    #[arg(short, long, value_name = "U8")]
    min_mapq: Option<MappingQuality>,

    /// Do not count supplementary alignments.
    #[arg(long)]
    no_supplementary: bool,

    /// Do count secondary alignments.
    #[arg(long)]
    count_secondary: bool,

    /// Do count duplicates.
    #[arg(long)]
    count_duplicates: bool,
}

/// Main function for the `ngs derive junction-annotation` subcommand.
pub fn derive(args: JunctionAnnotationArgs) -> anyhow::Result<()> {
    info!("Starting derive junction-annotation subcommand.");

    let mut exons = compute::ExonSets {
        starts: HashMap::new(),
        ends: HashMap::new(),
    };

    // (1) Parse the GFF file and collect all exon features.
    debug!("Reading all records in GFF.");
    let mut gff = formats::gff::open(&args.features_gff)
        .with_context(|| format!("opening GFF file: {}", args.features_gff.display()))?;

    let mut exon_records = Vec::new();
    for result in gff.records() {
        let record = result.unwrap();
        if record.ty() != args.exon_feature_name {
            continue;
        }
        exon_records.push(record);
    }
    debug!("Read {} exon records.", exon_records.len());

    debug!("Processing GFF exon features.");
    for record in &exon_records {
        let seq_name = record.reference_sequence_name();
        let start = record.start();
        let end = record.end().checked_add(1).unwrap(); // TODO: why +1? It works.

        exons.starts.entry(seq_name).or_default().insert(start);
        exons.ends.entry(seq_name).or_default().insert(end);
    }

    debug!("Done reading GFF.");

    // (1.5) Initialize variables (including opening the BAM).
    let mut counter = RecordCounter::default();
    let mut results = JunctionAnnotationResults::default();
    let params = compute::JunctionAnnotationParameters {
        min_intron_length: args.min_intron_length,
        min_read_support: args.min_read_support,
        min_mapq: args.min_mapq,
        no_supplementary: args.no_supplementary,
        count_secondary: args.count_secondary,
        count_duplicates: args.count_duplicates,
    };

    let ParsedBAMFile {
        mut reader, header, ..
    } = formats::bam::open_and_parse(args.src, IndexCheck::None)?;

    // (2) Process each record in the BAM file.
    for result in reader.records(&header.parsed) {
        let record = result?;
        compute::process(&record, &exons, &header.parsed, &params, &mut results)?;
        counter.inc();
    }

    info!(
        "Processed {} records.",
        counter.get().to_formatted_string(&Locale::en)
    );

    // (3) Summarize found junctions.
    compute::summarize(&mut results, &params);

    // (4) Print the output to stdout as JSON (more support for different output
    // types may be added in the future, but for now, only JSON).
    let output = serde_json::to_string_pretty(&results).unwrap();
    println!("{}", output);

    anyhow::Ok(())
}
