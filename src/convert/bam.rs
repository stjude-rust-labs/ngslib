//! Conversions from a BAM file to other next-generation sequencing file formats.

use std::io;
use std::num::NonZeroUsize;
use std::path::PathBuf;

use anyhow::Context;
use noodles::cram;
use noodles::fasta;
use noodles::sam;
use noodles::sam::alignment::Record;
use noodles::sam::header::record::value::map::ReferenceSequence;
use noodles::sam::header::record::value::Map;
use tokio::fs::File;
use tracing::info;

use crate::utils::args::NumberOfRecords;
use crate::utils::display::RecordCounter;
use crate::utils::formats;
use crate::utils::formats::bam::ParsedAsyncBAMFile;
use crate::utils::formats::utils::IndexCheck;

/// Converts a BAM file to a SAM file in an asyncronous fashion.
pub async fn to_sam_async(
    from: PathBuf,
    to: PathBuf,
    max_records: NumberOfRecords,
) -> anyhow::Result<()> {
    // (1) Open and parse the BAM file.
    let ParsedAsyncBAMFile {
        mut reader, header, ..
    } = formats::bam::open_and_parse_async(from, IndexCheck::None)
        .await
        .with_context(|| "opening BAM input file")?;

    // (2) Open the SAM file writer.
    let handle = File::create(to)
        .await
        .with_context(|| "creating SAM output file")?;
    let mut writer = sam::AsyncWriter::new(handle);

    // (3) Write the header.
    writer
        .write_header(&header.parsed)
        .await
        .with_context(|| "writing SAM header")?;

    let mut counter = RecordCounter::default();
    let mut record = Record::default();

    // (4) Write each record in the BAM file to the SAM file.
    while reader
        .read_record(&header.parsed, &mut record)
        .await
        .with_context(|| "reading BAM record")?
        != 0
    {
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

/// Converts a BAM file to a CRAM file in an asyncronous fashion.
pub async fn to_cram_async(
    from: PathBuf,
    to: PathBuf,
    fasta: PathBuf,
    max_records: NumberOfRecords,
) -> anyhow::Result<()> {
    // (1) Open and parse the BAM file.
    let ParsedAsyncBAMFile {
        mut reader,
        mut header,
        ..
    } = formats::bam::open_and_parse_async(from, IndexCheck::None)
        .await
        .with_context(|| "opening BAM input file")?;

    // (2) Builds the FASTA repository and associated index.
    let mut fasta_reader = formats::fasta::open(fasta)?;
    let records: Vec<fasta::Record> = fasta_reader
        .records()
        .collect::<io::Result<Vec<fasta::Record>>>()
        .with_context(|| "opening reference FASTA and associated index")?;

    // (3) Modifies the existing BAM header to include the reference sequences provided.
    let reference_sequences = header.parsed.reference_sequences_mut();
    for record in records.iter() {
        let name_as_string = record.name().to_owned();
        let name = name_as_string
            .parse()
            .with_context(|| "parsing reference FASTA sequence name")?;
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

    // (6) Write each record in the BAM file to the CRAM file.
    info!("Writing records to CRAM file.");
    while reader
        .read_record(&header.parsed, &mut record)
        .await
        .with_context(|| "reading BAM record")?
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
