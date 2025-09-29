//! Cell barcode processing functionality

use anyhow::Result;
use flate2::read::GzDecoder;
use rustc_hash::FxHashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::Arc;

/// Processor for cell barcodes
#[derive(Debug, Clone)]
pub struct BarcodeProcessor {
    /// Set of valid barcodes
    valid_barcodes: Arc<FxHashSet<String>>,
}

impl BarcodeProcessor {
    /// Create a new BarcodeProcessor from a file containing barcodes
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)?;
        let reader: Box<dyn BufRead> = if path
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("gz"))
            .unwrap_or(false)
        {
            Box::new(BufReader::with_capacity(256 * 1024, GzDecoder::new(file)))
        } else {
            Box::new(BufReader::with_capacity(256 * 1024, file))
        };

        Self::from_reader(reader)
    }

    /// Check if a barcode is valid
    #[inline]
    pub fn is_valid(&self, barcode: &str) -> bool {
        // Since barcodes are cached in an Arc<HashSet>, this is already efficient
        self.valid_barcodes.contains(barcode)
    }

    /// Get the number of valid barcodes
    pub fn len(&self) -> usize {
        self.valid_barcodes.len()
    }

    /// Check if the barcode set is empty
    pub fn is_empty(&self) -> bool {
        self.valid_barcodes.is_empty()
    }

    /// Get all valid barcodes as a sorted vector
    pub fn barcodes(&self) -> Vec<String> {
        let mut barcodes: Vec<_> = self.valid_barcodes.iter().cloned().collect();
        barcodes.sort();
        barcodes
    }

    fn from_reader(reader: Box<dyn BufRead>) -> Result<Self> {
        let mut valid_barcodes = FxHashSet::default();
        let mut line = String::with_capacity(64);

        for result in reader.lines() {
            line.clear();
            match result {
                Ok(line_content) => {
                    let barcode = line_content.trim();
                    if !barcode.is_empty() {
                        let clean_barcode = barcode.split('-').next().unwrap_or(barcode);
                        valid_barcodes.insert(clean_barcode.to_string());
                    }
                }
                Err(e) => return Err(e.into()),
            }
        }

        valid_barcodes.shrink_to_fit();
        Ok(BarcodeProcessor {
            valid_barcodes: Arc::new(valid_barcodes),
        })
    }
}
