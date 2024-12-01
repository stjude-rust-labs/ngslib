//! Conversions from a SAM file to other next-generation sequencing file formats.

use std::io;
use std::num::NonZeroUsize;
use std::path::PathBuf;

use anyhow::Context;
use noodles::bam;
use noodles::bgzf;
use noodles::bgzf::writer::CompressionLevel;
use noodles::cram;
use noodles::fasta;
use noodles::sam::alignment::Record;
use noodles::sam::header::record::value::map::ReferenceSequence;
use noodles::sam::header::record::value::Map;
use tokio::fs::File;
use tracing::info;

use crate::utils::args::CompressionStrategy;
use crate::utils::args::NumberOfRecords;
use crate::utils::display::RecordCounter;
use crate::utils::formats;
use crate::utils::formats::sam::ParsedAsyncSAMFile;

/// Converts a SAM file to a BAM file in an asyncronous fashion.
pub async fn to_bam_async(
    from: PathBuf,
    to: PathBuf,
    max_records: NumberOfRecords,
    compression_strategy: CompressionStrategy,
) -> anyhow::Result<()> {
    // (1) Open and parse the BAM file.
    let ParsedAsyncSAMFile { mut reader, header } = formats::sam::open_and_parse_async(from)
        .await
        .with_context(|| "opening SAM input file")?;

    // (2) Determine the compression level.
    let compression_level: CompressionLevel = compression_strategy.into();

    // (3) Open the SAM file writer.
    let mut writer = File::create(to)
        .await
        .map(|f| {
            bgzf::r#async::writer::Builder::default()
                .set_compression_level(compression_level)
                .build_with_writer(f)
        })
        .map(bam::AsyncWriter::from)
        .with_context(|| "opening BAM output file")?;

    // (4) Write the header and the reference sequences.
    writer
        .write_header(&header.parsed)
        .await
        .with_context(|| "writing BAM header")?;
    writer
        .write_reference_sequences(header.parsed.reference_sequences())
        .await
        .with_context(|| "writing BAM reference sequences")?;

    let mut counter = RecordCounter::default();
    let mut record = Record::default();

    // (5) Write each record in the BAM file to the SAM file.
    while reader
        .read_record(&header.parsed, &mut record)
        .await
        .with_context(|| "reading SAM record")?
        != 0
    {
        writer
            .write_alignment_record(&header.parsed, &record)
            .await
            .with_context(|| "writing BAM record")?;

        counter.inc();

        if counter.time_to_break(&max_records) {
            break;
        }
    }

    // (6) Shutdown the async writer.
    writer
        .shutdown()
        .await
        .with_context(|| "shutting down BAM writer")?;

    Ok(())
}

/// Converts a SAM file to a CRAM file in an asyncronous fashion.
pub async fn to_cram_async(
    from: PathBuf,
    to: PathBuf,
    fasta: PathBuf,
    max_records: NumberOfRecords,
) -> anyhow::Result<()> {
    // (1) Open and parse the SAM file.
    let ParsedAsyncSAMFile {
        mut reader,
        mut header,
        ..
    } = formats::sam::open_and_parse_async(from)
        .await
        .with_context(|| "reading SAM input file")?;

    // (2) Builds the FASTA repository and associated index.
    let mut fasta_reader =
        formats::fasta::open(fasta).with_context(|| "opening reference FASTA file")?;
    let records: Vec<fasta::Record> = fasta_reader
        .records()
        .collect::<io::Result<Vec<fasta::Record>>>()
        .with_context(|| "reading reference FASTA records")?;

    // (3) Modifies the existing SAM header to include the reference sequences provided.
    let reference_sequences = header.parsed.reference_sequences_mut();
    for record in records.iter() {
        let name_as_string = record.name().to_owned();
        let name = name_as_string
            .parse()
            .with_context(|| "parsing reference sequence name")?;
        let length = record.sequence().len();

        let reference_sequence = Map::<ReferenceSequence>::new(
            NonZeroUsize::try_from(length)
                .with_context(|| "converting reference sequence length to NonZeroUsize")?,
        );
        reference_sequences.insert(name, reference_sequence);
    }

    let repository = fasta::Repository::new(records);

    // (4) Open the CRAM file writer.
    let handle = File::create(to)
        .await
        .with_context(|| "opening CRAM output file")?;
    let mut writer = cram::r#async::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_with_writer(handle);

    // (5) Write the CRAM file.
    info!("Writing the file definition and header to CRAM file.");

    writer
        .write_file_definition()
        .await
        .with_context(|| "writing CRAM file definition")?;
    writer
        .write_file_header(&header.parsed)
        .await
        .with_context(|| "writing CRAM file header")?;

    let mut counter = RecordCounter::default();
    let mut record = Record::default();

    // (6) Write each record in the SAM file to the CRAM file.
    info!("Writing records to CRAM file.");
    while reader
        .read_record(&header.parsed, &mut record)
        .await
        .with_context(|| "reading SAM record")?
        != 0
    {
        let cram_record = cram::Record::try_from_alignment_record(&header.parsed, &record)?;
        writer
            .write_record(&header.parsed, cram_record)
            .await
            .with_context(|| "writing CRAM record")?;

        counter.inc();

        if counter.time_to_break(&max_records) {
            break;
        }
    }

    writer
        .shutdown(&header.parsed)
        .await
        .with_context(|| "shutting down CRAM writer")?;
    Ok(())
}
