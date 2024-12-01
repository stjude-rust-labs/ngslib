//! BAM indexing
//!
//! NOTICE: this was taken almost verbatim from @zaeleus's excellent example in
//! noodles. You can find that source code in the `bam_index.rs` example of the
//! `noodles_bam` crate at commit 0a709087934d.

use anyhow::bail;
use anyhow::Context;
use noodles::bam::bai;
use noodles::csi::index::reference_sequence::bin::Chunk;
use noodles::sam;
use noodles::sam::alignment::Record;
use noodles::sam::header::record::value::map::header::SortOrder;
use std::fs::File;
use std::path::PathBuf;
use tracing::debug;
use tracing::info;

use crate::utils::display::RecordCounter;
use crate::utils::formats::bam::ParsedBAMFile;
use crate::utils::formats::utils::IndexCheck;

//==================================//
// Individual indexing methods: BAM //
//==================================//

/// Checks if a SAM header indicates the file is coordinate sorted.
fn is_coordinate_sorted(header: &sam::Header) -> bool {
    if let Some(hdr) = header.header() {
        if let Some(sort_order) = hdr.sort_order() {
            return sort_order == SortOrder::Coordinate;
        }
    }

    false
}

/// Main method for BAM indexing.
pub fn index(src: PathBuf) -> anyhow::Result<()> {
    info!("indexing BAM file");

    // (1) Open and parse the BAM file.
    let ParsedBAMFile {
        mut reader,
        header,
        index_path,
        ..
    } = crate::utils::formats::bam::open_and_parse(src, IndexCheck::None)?;

    // (2) Calculate where the BAM index should go and check if a file is
    // already there. Error out if so.
    if index_path.exists() {
        bail!(
            "refusing to overwrite existing index file: {}. Please delete \
                and rerun if you'd like to replace it.",
            index_path.display()
        );
    }

    // (3) Check if the BAM is coordinate sorted. If it isn't, then return an error.
    debug!("checking the BAM is coordinate sorted");
    if !is_coordinate_sorted(&header.parsed) {
        bail!("the input BAM must be coordinate-sorted to be indexed");
    }

    // (5) Build the BAM index.
    debug!("building the BAM index");
    let mut record = Record::default();
    let mut builder = bai::Index::builder();
    let mut start_position = reader.virtual_position();

    let mut counter = RecordCounter::default();

    loop {
        match reader.read_record(&header.parsed, &mut record) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => bail!("failed to read record: {}", e),
        }

        let end_position = reader.virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        builder
            .add_record(&record, chunk)
            .with_context(|| "building BAM index")?;

        start_position = end_position;
        counter.inc();
    }

    debug!("building BAM index");
    let index = builder.build(header.parsed.reference_sequences().len());

    // (6) Write the index to disk.
    debug!("writing the BAM index to disk");
    let mut writer = File::create(index_path)
        .map(bai::Writer::new)
        .with_context(|| "creating BAM index output file")?;

    writer
        .write_header()
        .with_context(|| "writing BAM header")?;
    writer
        .write_index(&index)
        .with_context(|| "writing BAM index")?;

    Ok(())
}
