//! RNA editing type definitions with optimized strand-aware logic
//!
//! This module provides the core definitions and functionality for handling
//! different types of RNA editing events in a strand-aware manner. The
//! strand-aware logic ensures proper processing of editing events on both
//! positive and negative DNA strands by considering the complementary base
//! pairing rules.

use crate::core::error::{RedicatError, Result};
use log::{info, warn};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str::FromStr;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EditingType {
    AG,
    AC,
    AT,
    CA,
    CG,
    CT,
}

impl FromStr for EditingType {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "ag" => Ok(EditingType::AG),
            "ac" => Ok(EditingType::AC),
            "at" => Ok(EditingType::AT),
            "ca" => Ok(EditingType::CA),
            "cg" => Ok(EditingType::CG),
            "ct" => Ok(EditingType::CT),
            _ => Err(format!(
                "Invalid editing type: {}. Valid types: ag, ac, at, ca, cg, ct",
                s
            )),
        }
    }
}

impl fmt::Display for EditingType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            EditingType::AG => write!(f, "ag"),
            EditingType::AC => write!(f, "ac"),
            EditingType::AT => write!(f, "at"),
            EditingType::CA => write!(f, "ca"),
            EditingType::CG => write!(f, "cg"),
            EditingType::CT => write!(f, "ct"),
        }
    }
}

impl EditingType {
    /// Get the reference and alternate bases for this editing type
    pub fn to_bases(&self) -> (char, char) {
        match self {
            EditingType::AG => ('A', 'G'),
            EditingType::AC => ('A', 'C'),
            EditingType::AT => ('A', 'T'),
            EditingType::CA => ('C', 'A'),
            EditingType::CG => ('C', 'G'),
            EditingType::CT => ('C', 'T'),
        }
    }

    /// Get allowed reference bases for both strands (strand-aware)
    pub fn get_strand_aware_ref_bases(&self) -> [char; 2] {
        match self {
            EditingType::AG => ['A', 'T'],
            EditingType::AC => ['A', 'T'],
            EditingType::AT => ['A', 'T'],
            EditingType::CA => ['C', 'G'],
            EditingType::CG => ['C', 'G'],
            EditingType::CT => ['C', 'G'],
        }
    }

    /// Get the corresponding alt base for a given ref base considering strand
    pub fn get_alt_base_for_ref(&self, ref_base: char) -> char {
        match self {
            EditingType::AG => match ref_base {
                'A' => 'G',
                'T' => 'C',
                _ => 'N',
            },
            EditingType::AC => match ref_base {
                'A' => 'C',
                'T' => 'G',
                _ => 'N',
            },
            EditingType::AT => match ref_base {
                'A' => 'T',
                'T' => 'A',
                _ => 'N',
            },
            EditingType::CA => match ref_base {
                'C' => 'A',
                'G' => 'T',
                _ => 'N',
            },
            EditingType::CG => match ref_base {
                'C' => 'G',
                'G' => 'C',
                _ => 'N',
            },
            EditingType::CT => match ref_base {
                'C' => 'T',
                'G' => 'A',
                _ => 'N',
            },
        }
    }

    /// Check if a reference base is valid for this editing type
    pub fn is_valid_ref_base(&self, base: char) -> bool {
        self.get_strand_aware_ref_bases().contains(&base)
    }

    /// Get all valid editing types
    pub fn all_types() -> Vec<EditingType> {
        vec![
            EditingType::AG,
            EditingType::AC,
            EditingType::AT,
            EditingType::CA,
            EditingType::CG,
            EditingType::CT,
        ]
    }
}

/// Optimized REDIPortal loading with parallel processing and better memory management
pub fn load_rediportal_parallel(path: &str) -> Result<HashMap<String, u8>> {
    info!("Loading REDIPortal from: {}", path);

    if !std::path::Path::new(path).exists() {
        return Err(RedicatError::FileNotFound(format!(
            "REDIPortal file not found: {}",
            path
        )));
    }

    let file = File::open(path).map_err(|e| {
        RedicatError::FileNotFound(format!("Failed to open REDIPortal file {}: {}", path, e))
    })?;

    let reader: Box<dyn BufRead + Send> = if path.ends_with(".gz") {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let lines: Vec<String> = reader
        .lines()
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(RedicatError::Io)?;

    if lines.is_empty() {
        return Err(RedicatError::EmptyData(
            "REDIPortal file is empty".to_string(),
        ));
    }

    info!("Read {} lines from REDIPortal file", lines.len());

    let editing_sites: HashMap<String, u8> = lines
        .par_iter()
        .skip(1)
        .filter_map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();

            if fields.len() >= 2 {
                match fields[1].parse::<u64>() {
                    Ok(_pos) => {
                        let key = format!("{}:{}", fields[0], fields[1]);
                        Some((key, 1u8))
                    }
                    Err(_) => {
                        warn!("Invalid position in line: {}", line);
                        None
                    }
                }
            } else {
                warn!("Invalid line format (insufficient columns): {}", line);
                None
            }
        })
        .collect();

    if editing_sites.is_empty() {
        return Err(RedicatError::EmptyData(
            "No valid editing sites found in REDIPortal file".to_string(),
        ));
    }

    info!(
        "Loaded {} editing sites from REDIPortal",
        editing_sites.len()
    );

    if log::log_enabled!(log::Level::Debug) {
        let sample_sites: Vec<&String> = editing_sites.keys().take(5).collect();
        log::debug!("Sample editing sites: {:?}", sample_sites);
    }

    Ok(editing_sites)
}

/// Load REDIPortal with additional validation and filtering options
pub fn load_rediportal_with_filters(
    path: &str,
    min_chromosome_length: usize,
    allowed_chromosomes: Option<&[&str]>,
) -> Result<HashMap<String, u8>> {
    info!("Loading REDIPortal with filters from: {}", path);

    let mut editing_sites = load_rediportal_parallel(path)?;
    let original_count = editing_sites.len();

    if let Some(allowed_chroms) = allowed_chromosomes {
        editing_sites.retain(|key, _| {
            if let Some(chr) = key.split(':').next() {
                allowed_chroms.contains(&chr)
            } else {
                false
            }
        });
        info!(
            "Filtered by allowed chromosomes: {} -> {} sites",
            original_count,
            editing_sites.len()
        );
    }

    if min_chromosome_length > 0 {
        let before_filter = editing_sites.len();
        editing_sites.retain(|key, _| {
            if let Some(chr) = key.split(':').next() {
                chr.len() >= min_chromosome_length
            } else {
                false
            }
        });
        info!(
            "Filtered by chromosome name length (>= {}): {} -> {} sites",
            min_chromosome_length,
            before_filter,
            editing_sites.len()
        );
    }

    if editing_sites.is_empty() {
        return Err(RedicatError::EmptyData(
            "No editing sites remain after filtering".to_string(),
        ));
    }

    info!("Final filtered editing sites: {}", editing_sites.len());
    Ok(editing_sites)
}

/// Validate REDIPortal file format without loading all data
pub fn validate_rediportal_format(path: &str) -> Result<()> {
    info!("Validating REDIPortal file format: {}", path);

    let file = File::open(path).map_err(|e| {
        RedicatError::FileNotFound(format!("Failed to open REDIPortal file {}: {}", path, e))
    })?;

    let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut lines = reader.lines();

    let header = lines
        .next()
        .ok_or_else(|| RedicatError::EmptyData("REDIPortal file is empty".to_string()))??;

    info!("Header: {}", header);

    let mut valid_lines = 0;
    let mut invalid_lines = 0;
    const MAX_CHECK_LINES: usize = 100;

    for (i, line_result) in lines.enumerate().take(MAX_CHECK_LINES) {
        match line_result {
            Ok(line) => {
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() >= 2 && fields[1].parse::<u64>().is_ok() {
                    valid_lines += 1;
                } else {
                    invalid_lines += 1;
                    if invalid_lines <= 5 {
                        warn!("Invalid line {}: {}", i + 2, line);
                    }
                }
            }
            Err(e) => {
                return Err(RedicatError::Io(e));
            }
        }
    }

    if valid_lines == 0 {
        return Err(RedicatError::InvalidInput(
            "No valid data lines found in REDIPortal file".to_string(),
        ));
    }

    info!(
        "Validation complete: {} valid lines, {} invalid lines (checked first {} lines)",
        valid_lines, invalid_lines, MAX_CHECK_LINES
    );

    if invalid_lines > valid_lines / 2 {
        warn!("High proportion of invalid lines detected. Please check file format.");
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_editing_type_parsing() {
        assert_eq!(EditingType::from_str("ag").unwrap(), EditingType::AG);
        assert_eq!(EditingType::from_str("AG").unwrap(), EditingType::AG);
        assert_eq!(EditingType::from_str("Ag").unwrap(), EditingType::AG);
        assert!(EditingType::from_str("invalid").is_err());
    }

    #[test]
    fn test_strand_aware_bases() {
        let ag_type = EditingType::AG;
        assert_eq!(ag_type.get_strand_aware_ref_bases(), ['A', 'T']);
        assert_eq!(ag_type.get_alt_base_for_ref('A'), 'G');
        assert_eq!(ag_type.get_alt_base_for_ref('T'), 'C');
        assert_eq!(ag_type.get_alt_base_for_ref('G'), 'N');
    }

    #[test]
    fn test_editing_type_display() {
        assert_eq!(EditingType::AG.to_string(), "ag");
        assert_eq!(EditingType::CT.to_string(), "ct");
    }

    #[test]
    fn test_valid_ref_base() {
        let ag_type = EditingType::AG;
        assert!(ag_type.is_valid_ref_base('A'));
        assert!(ag_type.is_valid_ref_base('T'));
        assert!(!ag_type.is_valid_ref_base('G'));
        assert!(!ag_type.is_valid_ref_base('C'));
    }

    #[test]
    fn test_all_types_contains_expected_entries() {
        let types = EditingType::all_types();
        assert!(types.contains(&EditingType::AG));
        assert!(types.contains(&EditingType::CT));
        assert_eq!(types.len(), 6);
    }
}
