//! Cell barcode processing functionality

use anyhow::Result;
use flate2::read::GzDecoder;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::Arc;

/// Processor for cell barcodes
#[derive(Debug, Clone)]
pub struct BarcodeProcessor {
    ordered_barcodes: Arc<Vec<String>>,
    barcode_to_id: Arc<FxHashMap<String, u32>>,
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
        self.barcode_to_id.contains_key(barcode)
    }

    /// Lookup the numeric identifier for a barcode if present.
    #[inline]
    pub fn id_of(&self, barcode: &str) -> Option<u32> {
        self.barcode_to_id.get(barcode).copied()
    }

    /// Retrieve a barcode string by numeric identifier.
    #[inline]
    pub fn barcode_by_id(&self, id: u32) -> Option<&str> {
        self.ordered_barcodes.get(id as usize).map(|s| s.as_str())
    }

    /// Get the number of valid barcodes
    pub fn len(&self) -> usize {
        self.ordered_barcodes.len()
    }

    /// Check if the barcode set is empty
    pub fn is_empty(&self) -> bool {
        self.ordered_barcodes.is_empty()
    }

    /// Get all valid barcodes (preserving whitelist order)
    pub fn barcodes(&self) -> Vec<String> {
        self.ordered_barcodes.as_ref().clone()
    }

    /// Borrow the ordered whitelist without cloning.
    pub fn ordered_barcodes(&self) -> &[String] {
        self.ordered_barcodes.as_ref()
    }

    fn from_reader(reader: Box<dyn BufRead>) -> Result<Self> {
        let mut ordered = Vec::with_capacity(1024);
        let mut index = FxHashMap::default();
        let mut line = String::with_capacity(64);

        for result in reader.lines() {
            line.clear();
            match result {
                Ok(line_content) => {
                    let barcode = line_content.trim();
                    if !barcode.is_empty() {
                        let clean_barcode = barcode.split('-').next().unwrap_or(barcode);
                        if !index.contains_key(clean_barcode) {
                            let id = ordered.len() as u32;
                            ordered.push(clean_barcode.to_string());
                            index.insert(clean_barcode.to_string(), id);
                        }
                    }
                }
                Err(e) => return Err(e.into()),
            }
        }

        ordered.shrink_to_fit();
        index.shrink_to_fit();
        Ok(BarcodeProcessor {
            ordered_barcodes: Arc::new(ordered),
            barcode_to_id: Arc::new(index),
        })
    }

    /// Construct a processor from an explicit list of barcodes, preserving order.
    pub fn from_vec(barcodes: Vec<String>) -> Self {
        let mut index = FxHashMap::with_capacity_and_hasher(barcodes.len(), Default::default());
        for (i, barcode) in barcodes.iter().enumerate() {
            index.insert(barcode.clone(), i as u32);
        }

        BarcodeProcessor {
            ordered_barcodes: Arc::new(barcodes),
            barcode_to_id: Arc::new(index),
        }
    }
}
