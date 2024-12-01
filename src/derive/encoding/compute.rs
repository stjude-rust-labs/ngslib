//! Module holding the logic for computing the quality score encoding.

use anyhow::bail;
use serde::Serialize;
use std::collections::HashSet;

const MAX_VALID_PHRED_SCORE: u8 = 93;
const SANGER_MIN: u8 = 0;
const ILLUMINA_1_0_MIN: u8 = 26;
const ILLUMINA_1_3_MIN: u8 = 31;

/// Struct holding the final results for an `ngs derive encoding` subcommand
/// call.
#[derive(Debug, Serialize)]
pub struct DerivedEncodingResult {
    /// Whether or not the `ngs derive encoding` subcommand succeeded.
    pub succeeded: bool,

    /// The detected quality score encoding, if derivable.
    pub encoding: Option<String>,

    /// The minimum quality score observed.
    pub observed_min: u8,

    /// The maximum quality score observed.
    pub observed_max: u8,
}

impl DerivedEncodingResult {
    /// Creates a new [`DerivedEncodingResult`].
    pub fn new(
        succeeded: bool,
        encoding: Option<String>,
        observed_min: u8,
        observed_max: u8,
    ) -> Self {
        DerivedEncodingResult {
            succeeded,
            encoding,
            observed_min,
            observed_max,
        }
    }
}

/// Main method to evaluate the observed quality scores and
/// return a result for the derived encoding. This may fail, and the
/// resulting [`DerivedEncodingResult`] should be evaluated accordingly.
pub fn predict(score_set: HashSet<u8>) -> Result<DerivedEncodingResult, anyhow::Error> {
    if score_set.is_empty() {
        bail!("No quality scores were detected in the file.");
    }

    let observed_min = *score_set.iter().min().unwrap();
    let observed_max = *score_set.iter().max().unwrap();

    let mut result = DerivedEncodingResult::new(false, None, observed_min, observed_max);

    if observed_max > MAX_VALID_PHRED_SCORE {
        return anyhow::Ok(result);
    }
    match observed_min {
        ILLUMINA_1_3_MIN..=MAX_VALID_PHRED_SCORE => {
            result.succeeded = true;
            result.encoding = Some("Illumina 1.3".to_string());
        }
        ILLUMINA_1_0_MIN..=MAX_VALID_PHRED_SCORE => {
            result.succeeded = true;
            result.encoding = Some("Illumina 1.0".to_string());
        }
        SANGER_MIN..=MAX_VALID_PHRED_SCORE => {
            result.succeeded = true;
            result.encoding = Some("Sanger/Illumina 1.8".to_string());
        }
        _ => unreachable!(),
    }

    anyhow::Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_predict_illumina_1_3() {
        let mut score_set: HashSet<u8> = HashSet::new();
        score_set.insert(40);
        score_set.insert(41);
        score_set.insert(42);
        score_set.insert(43);
        score_set.insert(44);
        score_set.insert(45);
        score_set.insert(46);
        score_set.insert(47);
        score_set.insert(48);
        score_set.insert(49);
        score_set.insert(50);
        score_set.insert(51);
        score_set.insert(52);
        score_set.insert(53);
        score_set.insert(54);
        score_set.insert(55);
        score_set.insert(56);
        score_set.insert(57);
        score_set.insert(58);
        score_set.insert(59);
        score_set.insert(60);
        score_set.insert(61);
        score_set.insert(62);
        score_set.insert(63);
        score_set.insert(64);
        score_set.insert(65);
        score_set.insert(66);
        score_set.insert(67);
        score_set.insert(68);
        score_set.insert(69);
        score_set.insert(70);
        score_set.insert(71);
        score_set.insert(72);
        score_set.insert(73);
        score_set.insert(74);
        score_set.insert(75);
        score_set.insert(76);
        score_set.insert(77);
        score_set.insert(78);
        score_set.insert(79);
        score_set.insert(80);
        score_set.insert(81);
        score_set.insert(82);
        score_set.insert(83);
        score_set.insert(84);
        score_set.insert(85);
        score_set.insert(86);
        score_set.insert(87);
        score_set.insert(88);
        score_set.insert(89);
        score_set.insert(90);
        score_set.insert(91);
        score_set.insert(92);
        score_set.insert(93);

        let result = predict(score_set).unwrap();
        assert!(result.succeeded);
        assert_eq!(result.encoding, Some("Illumina 1.3".to_string()));
        assert_eq!(result.observed_min, 40);
        assert_eq!(result.observed_max, 93);
    }

    #[test]
    fn test_predict_illumina_1_0() {
        let mut score_set: HashSet<u8> = HashSet::new();
        score_set.insert(26);
        score_set.insert(27);
        score_set.insert(28);
        score_set.insert(29);
        score_set.insert(30);
        score_set.insert(31);
        score_set.insert(32);
        score_set.insert(33);
        score_set.insert(34);
        score_set.insert(35);
        score_set.insert(36);
        score_set.insert(37);
        score_set.insert(38);
        score_set.insert(39);
        score_set.insert(40);
        score_set.insert(41);
        score_set.insert(42);
        score_set.insert(43);
        score_set.insert(44);
        score_set.insert(45);
        score_set.insert(46);
        score_set.insert(47);
        score_set.insert(48);
        score_set.insert(49);
        score_set.insert(50);
        score_set.insert(51);
        score_set.insert(52);
        score_set.insert(53);
        score_set.insert(54);
        score_set.insert(55);
        score_set.insert(56);
        score_set.insert(57);
        score_set.insert(58);
        score_set.insert(59);
        score_set.insert(60);
        score_set.insert(61);
        score_set.insert(62);
        score_set.insert(63);
        score_set.insert(64);
        score_set.insert(65);
        score_set.insert(66);
        score_set.insert(67);
        score_set.insert(68);
        score_set.insert(69);
        score_set.insert(70);
        score_set.insert(71);
        score_set.insert(72);
        score_set.insert(73);
        score_set.insert(74);
        score_set.insert(75);
        score_set.insert(76);
        score_set.insert(77);
        score_set.insert(78);
        score_set.insert(79);
        score_set.insert(80);
        score_set.insert(81);
        score_set.insert(82);
        score_set.insert(83);
        score_set.insert(84);
        score_set.insert(85);
        score_set.insert(86);
        score_set.insert(87);
        score_set.insert(88);
        score_set.insert(89);
        score_set.insert(90);
        score_set.insert(91);
        score_set.insert(92);
        score_set.insert(93);

        let result = predict(score_set).unwrap();
        assert!(result.succeeded);
        assert_eq!(result.encoding, Some("Illumina 1.0".to_string()));
        assert_eq!(result.observed_min, 26);
        assert_eq!(result.observed_max, 93);
    }

    #[test]
    fn test_predict_sanger() {
        let mut score_set: HashSet<u8> = HashSet::new();
        score_set.insert(0);
        score_set.insert(1);
        score_set.insert(2);
        score_set.insert(3);
        score_set.insert(4);
        score_set.insert(5);
        score_set.insert(6);
        score_set.insert(7);
        score_set.insert(8);
        score_set.insert(9);
        score_set.insert(10);
        score_set.insert(11);
        score_set.insert(12);
        score_set.insert(13);
        score_set.insert(14);
        score_set.insert(15);
        score_set.insert(16);
        score_set.insert(17);
        score_set.insert(18);
        score_set.insert(19);
        score_set.insert(20);
        score_set.insert(21);
        score_set.insert(22);
        score_set.insert(23);
        score_set.insert(24);
        score_set.insert(25);
        score_set.insert(26);
        score_set.insert(27);
        score_set.insert(28);
        score_set.insert(29);
        score_set.insert(30);
        score_set.insert(31);
        score_set.insert(32);
        score_set.insert(33);
        score_set.insert(34);
        score_set.insert(35);
        score_set.insert(36);
        score_set.insert(37);
        score_set.insert(38);
        score_set.insert(39);
        score_set.insert(40);
        score_set.insert(41);
        score_set.insert(42);
        score_set.insert(43);
        score_set.insert(44);
        score_set.insert(45);
        score_set.insert(46);
        score_set.insert(47);
        score_set.insert(48);
        score_set.insert(49);
        score_set.insert(50);
        score_set.insert(51);
        score_set.insert(52);
        score_set.insert(53);
        score_set.insert(54);
        score_set.insert(55);
        score_set.insert(56);
        score_set.insert(57);
        score_set.insert(58);
        score_set.insert(59);
        score_set.insert(60);
        score_set.insert(61);
        score_set.insert(62);
        score_set.insert(63);
        score_set.insert(64);
        score_set.insert(65);
        score_set.insert(66);
        score_set.insert(67);
        score_set.insert(68);

        let result = predict(score_set).unwrap();
        assert!(result.succeeded);
        assert_eq!(result.encoding, Some("Sanger/Illumina 1.8".to_string()));
        assert_eq!(result.observed_min, 0);
        assert_eq!(result.observed_max, 68);
    }

    #[test]
    fn test_predict_fail() {
        let score_set: HashSet<u8> = HashSet::new();
        let result = predict(score_set);
        assert!(result.is_err());
    }

    #[test]
    fn test_predict_too_high_max_score() {
        let mut score_set: HashSet<u8> = HashSet::new();
        score_set.insert(94);
        let result = predict(score_set).unwrap();
        assert!(!result.succeeded);
        assert_eq!(result.encoding, None);
        assert_eq!(result.observed_min, 94);
        assert_eq!(result.observed_max, 94);
    }
}
