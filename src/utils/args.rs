//! Utilities related to the parsing of arguments.

use std::fmt::Display;
use std::num::NonZeroUsize;

use noodles::bgzf::writer::CompressionLevel;
use tracing::debug;

//===================//
// Number of Records //
//===================//

/// Utility enum to designate whether we are reviewing all records in the file
/// or just some of them.
#[derive(Clone, Debug)]
pub enum NumberOfRecords {
    /// Designates that we should review _all_ of the records in the file.
    All,

    /// Designates that we should review _some_ of the records in the file. The
    /// exact count of records is stored in the `usize`.
    Some(NonZeroUsize),
}

impl std::default::Default for NumberOfRecords {
    fn default() -> Self {
        Self::All
    }
}

impl std::fmt::Display for NumberOfRecords {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NumberOfRecords::All => write!(f, "all"),
            NumberOfRecords::Some(value) => write!(f, "{value}"),
        }
    }
}

/// An error type for parsing the number of records.
#[derive(Debug)]
pub enum NumberOfRecordsError {
    /// The number of records is invalid.
    Invalid(String),
}

impl std::fmt::Display for NumberOfRecordsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NumberOfRecordsError::Invalid(value) => write!(f, "invalid number of reads: {value}"),
        }
    }
}

impl std::error::Error for NumberOfRecordsError {}

impl std::str::FromStr for NumberOfRecords {
    type Err = NumberOfRecordsError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "all" => Ok(NumberOfRecords::All),
            _ => s
                .parse::<usize>()
                .map_err(|_| {
                    NumberOfRecordsError::Invalid(String::from(
                        "must be a positive, non-zero integer or 'all'",
                    ))
                })
                .and_then(|num| {
                    NonZeroUsize::new(num)
                        .ok_or_else(|| {
                            NumberOfRecordsError::Invalid(String::from(
                                "integers must be positive and non-zero",
                            ))
                        })
                        .map(NumberOfRecords::Some)
                }),
        }
    }
}

impl From<Option<NonZeroUsize>> for NumberOfRecords {
    fn from(num_records: Option<NonZeroUsize>) -> Self {
        match num_records {
            Some(n) => {
                debug!("Reading a maximum of {} records.", n);
                NumberOfRecords::Some(n)
            }
            None => {
                debug!("Reading all available records.");
                NumberOfRecords::All
            }
        }
    }
}

//======================//
// Compression Strategy //
//======================//

/// An enum representing the compression strategy to follow.
#[derive(clap::ValueEnum, Clone, Debug)]
pub enum CompressionStrategy {
    /// Compress the file as much as possible (maximum gzip compression).
    Best,

    /// Balance the compression level and the speed of the compression process.
    Balanced,

    /// Compress the file as quickly as possible (minimum gzip compression without
    /// turning off compression altogether).
    Fastest,
}

impl Display for CompressionStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Balanced => write!(f, "balanced"),
            Self::Best => write!(f, "best"),
            Self::Fastest => write!(f, "fastest"),
        }
    }
}

impl From<CompressionStrategy> for CompressionLevel {
    fn from(strategy: CompressionStrategy) -> Self {
        match strategy {
            CompressionStrategy::Balanced => CompressionLevel::try_from(6)
                .expect("can create a compression level with value of 6"),
            CompressionStrategy::Best => CompressionLevel::best(),
            CompressionStrategy::Fastest => CompressionLevel::fast(),
        }
    }
}

//=============//
// Arg Parsers //
//=============//

/// Utility method to parse command line floats and ensure they are
/// within the range [MIN, MAX].
pub fn arg_in_range(arg: f64, range: std::ops::RangeInclusive<f64>) -> anyhow::Result<f64> {
    match range.contains(&arg) {
        true => Ok(arg),
        false => anyhow::bail!(
            "Value must be between {} and {}",
            range.start(),
            range.end()
        ),
    }
}
