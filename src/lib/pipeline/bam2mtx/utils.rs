//! Utility functions for bam2mtx processing

use std::path::Path;

use anyhow::Result;
use csv::ReaderBuilder;
use log::info;

/// Read TSV file with genomic positions
/// Supports both .tsv and .tsv.gz formats
pub fn read_positions_from_tsv<P: AsRef<Path>>(
    path: P,
    chrom_col: &str,
    pos_col: &str,
) -> Result<Vec<(String, u64)>> {
    use flate2::read::GzDecoder;
    use std::fs::File;
    use std::io::BufReader;

    // Create reader based on file extension
    let file = File::open(&path)?;
    let buf_reader: Box<dyn std::io::BufRead> =
        if path.as_ref().extension().map_or(false, |ext| ext == "gz") {
            // Handle .tsv.gz files
            let decoder = GzDecoder::new(file);
            Box::new(BufReader::new(decoder))
        } else {
            // Handle regular .tsv files
            Box::new(BufReader::new(file))
        };

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(buf_reader);

    let headers = reader.headers()?.clone();
    let chrom_idx = headers
        .iter()
        .position(|h| h == chrom_col)
        .ok_or_else(|| anyhow::anyhow!("Column '{}' not found", chrom_col))?;
    let pos_idx = headers
        .iter()
        .position(|h| h == pos_col)
        .ok_or_else(|| anyhow::anyhow!("Column '{}' not found", pos_col))?;

    let mut positions = Vec::new();

    for result in reader.records() {
        let record = result?;
        let chrom = record
            .get(chrom_idx)
            .ok_or_else(|| anyhow::anyhow!("Missing chromosome value"))?;
        let pos_str = record
            .get(pos_idx)
            .ok_or_else(|| anyhow::anyhow!("Missing position value"))?;

        let pos: u64 = pos_str
            .parse()
            .map_err(|_| anyhow::anyhow!("Invalid position: {}", pos_str))?;

        positions.push((chrom.to_string(), pos));
    }

    Ok(positions)
}

/// Chunk positions into smaller batches for parallel processing
pub fn chunk_positions(positions: &[(String, u64)], chunk_size: usize) -> Vec<Vec<(String, u64)>> {
    positions
        .chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect()
}

/// Simple progress bar for CLI applications
pub struct ProgressBar {
    /// Total number of items to process
    total: usize,
    /// Number of items processed so far
    current: usize,
    /// Last percentage that was printed
    last_percent: usize,
}

impl ProgressBar {
    /// Create a new ProgressBar with the specified total count
    pub fn new(total: usize) -> Self {
        Self {
            total,
            current: 0,
            last_percent: 0,
        }
    }

    /// Increment the progress by the specified amount and update the display
    pub fn increment(&mut self, by: usize) {
        self.current += by;
        let percent = (self.current * 100) / self.total;

        if percent != self.last_percent {
            self.last_percent = percent;
            info!("\rProgress: {}% ({}/{})", percent, self.current, self.total);
            if self.current == self.total {
                info!("");
            }
        }
    }
}

/// Merge multiple PositionData results
pub fn merge_position_data(
    results: Vec<Vec<crate::pipeline::bam2mtx::processor::PositionData>>,
) -> Vec<crate::pipeline::bam2mtx::processor::PositionData> {
    let mut merged = Vec::new();

    for batch in results {
        merged.extend(batch);
    }

    merged
}
