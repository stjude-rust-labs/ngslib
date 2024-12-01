//! Conversions from a CRAM file to other next-generation sequencing file formats.

use std::path::PathBuf;

use anyhow::Context;
use futures::TryStreamExt;
use noodles::bam;
use noodles::bgzf;
use noodles::bgzf::writer::CompressionLevel;
use noodles::fasta;
use noodles::sam;
use tokio::fs::File;

use crate::utils::args::CompressionStrategy;
use crate::utils::args::NumberOfRecords;
use crate::utils::display::RecordCounter;
use crate::utils::formats;
use crate::utils::formats::cram::ParsedAsyncCRAMFile;
use crate::utils::formats::utils::IndexCheck;

/// Converts a CRAM file to a SAM file in an asyncronous fashion.
pub async fn to_sam_async(
    from: PathBuf,
    to: PathBuf,
    fasta: PathBuf,
    max_records: NumberOfRecords,
) -> anyhow::Result<()> {
    // (1) Open and parse the CRAM file.
    let ParsedAsyncCRAMFile {
        mut reader, header, ..
    } = formats::cram::open_and_parse_async(from, IndexCheck::None)
        .await
        .with_context(|| "opening CRAM file")?;

    // (2) Builds the FASTA repository and associated index.
    let repository = fasta::indexed_reader::Builder::default()
        .build_from_path(fasta)
        .map(fasta::repository::adapters::IndexedReader::new)
        .map(fasta::Repository::new)
        .with_context(|| "opening reference FASTA and associated index")?;

    // (3) Open the SAM file writer.
    let handle = File::create(to)
        .await
        .with_context(|| "opening SAM output file")?;
    let mut writer = sam::AsyncWriter::new(handle);

    // (4) Write the header.
    writer
        .write_header(&header.parsed)
        .await
        .with_context(|| "writing SAM header")?;

    // (5) Write each record in the CRAM file to the SAM file.
    let mut counter = RecordCounter::default();
    let mut records = reader.records(&repository, &header.parsed);

    while let Some(record) = records
        .try_next()
        .await
        .with_context(|| "reading CRAM record")?
    {
        let record = record
            .try_into_alignment_record(&header.parsed)
            .with_context(|| "parsing CRAM record")?;
        writer
            .write_alignment_record(&header.parsed, &record)
            .await
            .with_context(|| "writing SAM record")?;

        counter.inc();

        if counter.time_to_break(&max_records) {
            break;
        }
    }

    Ok(())
}

/// Converts a CRAM file to a BAM file in an asyncronous fashion.
pub async fn to_bam_async(
    from: PathBuf,
    to: PathBuf,
    fasta: PathBuf,
    max_records: NumberOfRecords,
    compression_strategy: CompressionStrategy,
) -> anyhow::Result<()> {
    // (1) Open and parse the CRAM file.
    let ParsedAsyncCRAMFile {
        mut reader, header, ..
    } = formats::cram::open_and_parse_async(from, IndexCheck::None)
        .await
        .with_context(|| "opening CRAM input file")?;

    // (2) Builds the FASTA repository and associated index.
    let repository = fasta::indexed_reader::Builder::default()
        .build_from_path(fasta)
        .map(fasta::repository::adapters::IndexedReader::new)
        .map(fasta::Repository::new)
        .with_context(|| "opening reference FASTA and associated index")?;

    // (3) Determine the compression level.
    let compression_level: CompressionLevel = compression_strategy.into();

    // (4) Open the BAM file writer.
    let mut writer = File::create(to)
        .await
        .map(|f| {
            bgzf::r#async::writer::Builder::default()
                .set_compression_level(compression_level)
                .build_with_writer(f)
        })
        .map(bam::AsyncWriter::from)
        .with_context(|| "opening BAM output file")?;

    // (5) Write the header and the reference sequences.
    writer
        .write_header(&header.parsed)
        .await
        .with_context(|| "writing BAM header")?;
    writer
        .write_reference_sequences(header.parsed.reference_sequences())
        .await
        .with_context(|| "writing BAM reference sequences")?;

    // (6) Write each record in the CRAM file to the BAM file.
    let mut counter = RecordCounter::default();
    let mut records = reader.records(&repository, &header.parsed);

    while let Some(record) = records
        .try_next()
        .await
        .with_context(|| "reading CRAM record")?
    {
        let record = record
            .try_into_alignment_record(&header.parsed)
            .with_context(|| "parsing CRAM record")?;
        writer
            .write_alignment_record(&header.parsed, &record)
            .await
            .with_context(|| "writing BAM record")?;

        counter.inc();

        if counter.time_to_break(&max_records) {
            break;
        }
    }

    // (7) Shutdown the async writer.
    writer
        .shutdown()
        .await
        .with_context(|| "shutting down BAM writer")?;

    Ok(())
}
