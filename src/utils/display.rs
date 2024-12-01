//! Utilities related to displaying things.

use std::fmt;

use num_format::Locale;
use num_format::ToFormattedString;
use tracing::info;

use crate::utils::args::NumberOfRecords;

/// Utility struct for displays percentages. The first item in the struct is the
/// numerator and the second item in the struct is the denominator.
pub struct PercentageFormat(pub u64, pub u64);

impl fmt::Display for PercentageFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.1 == 0 {
            f.write_str("N/A")
        } else {
            let (a, b) = (self.0 as f64, self.1 as f64);
            write!(f, "{:.2}%", a / b * 100.0)
        }
    }
}

/// Utility struct used to uniformly count and report the number of records processed.
pub struct RecordCounter {
    /// The number of records processed.
    count: usize,

    /// The number of records to log every.
    log_every: usize,
}

impl Default for RecordCounter {
    fn default() -> Self {
        RecordCounter {
            count: 0,
            log_every: 1_000_000,
        }
    }
}

impl RecordCounter {
    /// Creates a new `RecordCounter`.
    pub fn new(log_every: Option<usize>) -> Self {
        RecordCounter {
            count: 0,
            log_every: log_every.unwrap_or(1_000_000),
        }
    }

    /// Gets the current number of records counted via a copy.
    pub fn get(&self) -> usize {
        self.count
    }

    /// Increments the counter and reports the number of records processed (if
    /// appropriate).
    pub fn inc(&mut self) {
        self.count += 1;

        if self.count % self.log_every == 0 {
            info!(
                "  [*] Processed {} records.",
                self.count.to_formatted_string(&Locale::en),
            );
        }
    }

    /// A utility method that indicates whether a loop should break based on if the
    /// counter is greater than or equal to some limit. This is especially useful if you
    /// have an `Option<usize>` that indicates the maximum number of records to process
    /// (if it exists, otherwise it loops forever).
    pub fn time_to_break(&self, limit: &NumberOfRecords) -> bool {
        match limit {
            NumberOfRecords::Some(v) => self.count >= <usize>::from(*v),
            NumberOfRecords::All => false,
        }
    }
}
