//! Functionality related to the `ngs derive` subcommand itself.

pub mod encoding;
pub mod endedness;
pub mod instrument;
pub mod junction_annotation;
pub mod readlen;
pub mod strandedness;

use clap::Args;
use clap::Subcommand;

//===============//
// Command setup //
//===============//

/// Command line arguments for `ngs derive`.
#[derive(Args)]
pub struct DeriveArgs {
    /// The subcommand for `ngs derive`.
    #[command(subcommand)]
    pub subcommand: DeriveSubcommand,
}

/// All possible subcommands for `ngs derive`.
#[derive(Subcommand)]
pub enum DeriveSubcommand {
    /// Derives the quality score encoding used to produce the file.
    Encoding(self::encoding::DeriveEncodingArgs),

    /// Derives the endedness of the file.
    Endedness(self::endedness::DeriveEndednessArgs),

    /// Derives the instrument used to produce the file.
    Instrument(self::instrument::DeriveInstrumentArgs),

    /// Derives the read length of the file.
    Readlen(self::readlen::DeriveReadlenArgs),

    /// Derives the strandedness of the RNA-Seq file.
    Strandedness(self::strandedness::DeriveStrandednessArgs),

    /// Annotates junctions in the file.
    /// Note that, technically, this command doesn't derive anything—it will moved in the future to a better home.
    /// convenience.
    JunctionAnnotation(self::junction_annotation::JunctionAnnotationArgs),
}
