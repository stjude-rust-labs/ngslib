//! Functionality relating to the `ngs derive strandedness` subcommand itself.

use anyhow::{bail, Context};
use clap::Args;
use noodles::sam::record::MappingQuality;
use noodles::{bam, gff};
use rust_lapper::{Interval, Lapper};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::PathBuf;
use tracing::{debug, info};

use crate::derive::strandedness::compute::ParsedBAMFile;
use crate::derive::strandedness::{compute, results};
use crate::utils::formats;
use crate::utils::read_groups::validate_read_group_info;

/// Clap arguments for the `ngs derive strandedness` subcommand.
#[derive(Args)]
pub struct DeriveStrandednessArgs {
    /// Source BAM.
    #[arg(value_name = "BAM")]
    src: PathBuf,

    /// Features GFF file.
    #[arg(short = 'f', long, required = true, value_name = "PATH")]
    features_gff: PathBuf,

    /// When inconclusive, the test will repeat until this many tries have been reached.
    /// Evidence of previous attempts is saved and reused,
    /// leading to a larger sample size with multiple attempts.
    #[arg(long, value_name = "USIZE", default_value = "3")]
    max_tries: usize,

    /// Filter any genes that don't have at least `m` reads.
    #[arg(short = 'm', long, value_name = "USIZE", default_value = "10")]
    min_reads_per_gene: usize,

    /// How many genes to use as evidence in strandendess classification per try.
    /// This does not count genes which fail filtering
    /// due to `--min-reads-per-gene` or are discarded
    /// due to problematic Strand information in the GFF.
    /// Problematic Strand information is caused by contradictions between
    /// gene entries and overlapping exon entries.
    #[arg(short = 'n', long, value_name = "USIZE", default_value = "1000")]
    num_genes: usize,

    /// Minumum mapping quality for a record to be considered.
    /// Default behavior is to ignore MAPQ values,
    /// which allows reads with _missing_ MAPQs to be considered.
    /// Specify any u8 value (lower than 255) to enable this filter.
    /// Some aligners erroneously use 255 as the score for a uniquely mapped read;
    /// however, 255 is reserved by the spec for a missing MAPQ value.
    /// Therefore BAMs produced by aligners using 255 erroneously
    /// are not compatible with setting this option.
    #[arg(long, value_name = "U8")]
    min_mapq: Option<MappingQuality>,

    /// Consider all genes, not just protein coding genes.
    #[arg(long)]
    all_genes: bool,

    /// Name of the gene region feature for the gene model used.
    #[arg(long, value_name = "STRING", default_value = "gene")]
    gene_feature_name: String,

    /// Name of the exon region feature for the gene model used.
    #[arg(long, value_name = "STRING", default_value = "exon")]
    exon_feature_name: String,

    /// Do not count supplementary alignments.
    #[arg(long)]
    no_supplementary: bool,

    /// Do count secondary alignments.
    #[arg(long)]
    count_secondary: bool,

    /// Do count duplicates.
    #[arg(long)]
    count_duplicates: bool,

    /// Do count QC failed reads.
    #[arg(long)]
    count_qc_failed: bool,

    /// At most, evaluate this many genes
    /// per try. Default is 10 * --num-genes.
    #[arg(long, value_name = "USIZE")]
    max_genes_per_try: Option<usize>,
}

/// Main function for the `ngs derive strandedness` subcommand.
pub fn derive(args: DeriveStrandednessArgs) -> anyhow::Result<()> {
    info!("Starting derive strandedness subcommand.");

    // (1) Parse the GFF file and collect all gene and exon features.
    debug!("Reading all records in GFF.");
    let mut gff = formats::gff::open(&args.features_gff)
        .with_context(|| format!("opening GFF file: {}", args.features_gff.display()))?;

    let mut gene_records = Vec::new();
    let mut exon_records = Vec::new();
    let mut gene_metrics = results::GeneRecordMetrics::default();
    let mut exon_metrics = results::ExonRecordMetrics::default();
    for result in gff.records() {
        let record = result.unwrap();
        if record.ty() == args.gene_feature_name {
            gene_metrics.total += 1;

            // If --all-genes is set, don't check the gene type or biotype.
            // Otherwise, check the gene type or biotype and keep the record if it's protein coding.
            // If the record does not have a gene type or biotype, discard it.
            let mut keep_record = false;
            if !args.all_genes {
                let mut gene_type_value = None;
                for entry in record.attributes().as_ref() {
                    gene_type_value = match entry.key() {
                        "gene_type" => Some(entry.value()),    // Gencode
                        "gene_biotype" => Some(entry.value()), // ENSEMBL
                        "biotype" => Some(entry.value()),      // also ENSEMBL
                        _ => gene_type_value,
                    };
                }
                if let Some(gene_type_value) = gene_type_value {
                    if gene_type_value.to_lowercase().contains("protein") {
                        keep_record = true;
                        gene_metrics.protein_coding += 1;
                    }
                }
            }
            if !keep_record {
                continue;
            }

            // Make sure the gene record has a valid strand.
            let gene_strand = record.strand();
            if gene_strand != gff::record::Strand::Forward
                && gene_strand != gff::record::Strand::Reverse
            {
                gene_metrics.bad_strand += 1;
                continue;
            }

            gene_records.push(record);
        } else if record.ty() == args.exon_feature_name {
            exon_metrics.total += 1;
            exon_records.push(record);
        }
    }
    if gene_records.is_empty() {
        bail!("No gene records matched criteria. Check your GFF file and `--gene-feature-name` and `--all-genes` options.");
    }
    if exon_records.is_empty() {
        bail!("No exon records matched criteria. Check your GFF file and `--exon-feature-name` option.");
    }
    debug!(
        "Found {} gene records and {} exon records.",
        gene_records.len(),
        exon_records.len()
    );

    // (2) Parse exon features into proper data structure.
    debug!("Tabulating GFF exon features.");

    let mut exon_intervals: HashMap<&str, Vec<Interval<usize, compute::Strand>>> = HashMap::new();
    for record in &exon_records {
        let seq_name = record.reference_sequence_name();
        let start: usize = record.start().into();
        let stop: usize = record.end().into();
        let strand = record.strand();

        if strand != gff::record::Strand::Forward && strand != gff::record::Strand::Reverse {
            exon_metrics.bad_strand += 1;
            continue;
        }
        let strand = compute::Strand::try_from(strand).unwrap(); // above check guarantees safety

        exon_intervals.entry(seq_name).or_default().push(Interval {
            start,
            stop,
            val: strand,
        });
    }

    if exon_metrics.bad_strand == exon_metrics.total {
        bail!("All exons were discarded due to bad strand information. Check your GFF file.");
    }
    debug!(
        "{} exons were discarded due to bad strand information.",
        exon_metrics.bad_strand
    );

    let mut exons: HashMap<&str, Lapper<usize, compute::Strand>> = HashMap::new();
    for (seq_name, intervals) in exon_intervals {
        exons.insert(seq_name, Lapper::new(intervals));
    }

    debug!("Done reading GFF.");

    // (3) Initialize variables (including opening the BAM).
    let mut reader = File::open(&args.src)
        .map(bam::Reader::new)
        .with_context(|| format!("opening BAM file: {}", args.src.display()))?;
    let header = reader.read_header()?.parse()?;
    let index = bam::bai::read(args.src.with_extension("bam.bai")).with_context(|| {
        format!(
            "reading BAM index: {}",
            args.src.with_extension("bam.bai").display()
        )
    })?;

    let mut parsed_bam = ParsedBAMFile {
        reader,
        header,
        index,
    };

    let max_genes_per_try = args.max_genes_per_try.unwrap_or(args.num_genes * 10);

    let params = compute::StrandednessParams {
        num_genes: args.num_genes,
        max_genes_per_try,
        min_reads_per_gene: args.min_reads_per_gene,
        min_mapq: args.min_mapq,
        count_qc_failed: args.count_qc_failed,
        no_supplementary: args.no_supplementary,
        count_secondary: args.count_secondary,
        count_duplicates: args.count_duplicates,
    };

    let mut all_counts = compute::AllReadGroupsCounts {
        counts: HashMap::new(),
        found_rgs: HashSet::new(),
    };
    let mut metrics = results::RecordTracker {
        genes: gene_metrics,
        exons: exon_metrics,
        reads: results::ReadRecordMetrics::default(),
    };
    let mut result: Option<results::DerivedStrandednessResult> = None;

    // (4) Run the strandedness test.
    for try_num in 1..=args.max_tries {
        info!("Starting try {} of {}", try_num, args.max_tries);

        let attempt = compute::predict(
            &mut parsed_bam,
            &mut gene_records,
            &exons,
            &mut all_counts,
            &params,
            &mut metrics,
        )?;
        let success = attempt.succeeded;
        result = Some(attempt);
        if success {
            info!("Strandedness test succeeded.");
            break;
        } else {
            info!("Strandedness test inconclusive.");
        }
    }
    let mut result = result.unwrap();

    if !result.succeeded {
        info!(
            "Strandedness test still failed after {} tries.",
            args.max_tries
        );
    }

    let rgs_in_header_not_found =
        validate_read_group_info(&all_counts.found_rgs, &parsed_bam.header);
    let mut empty_rg_results = Vec::new();
    for rg in rgs_in_header_not_found {
        empty_rg_results.push(compute::predict_strandedness(
            &rg,
            &compute::Counts::default(),
        ));
    }
    result.read_groups.extend(empty_rg_results);

    // (5) Print the output to stdout as JSON (more support for different output
    // types may be added in the future, but for now, only JSON).
    let output = serde_json::to_string_pretty(&result).unwrap();
    println!("{}", output);

    anyhow::Ok(())
}
