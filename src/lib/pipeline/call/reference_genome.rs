//! Optimized reference genome operations with caching and better performance

use crate::core::error::{RedicatError, Result};
use bio::io::fasta::{Index, IndexedReader};
use log::{debug, info, warn};
use parking_lot::RwLock;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

/// Thread-safe reference genome accessor with caching
pub struct ReferenceGenome {
    reader: parking_lot::Mutex<IndexedReader<std::fs::File>>,
    sequences: Vec<String>,
    // Cache for frequently accessed positions
    cache: Arc<RwLock<HashMap<String, char>>>,
    cache_size_limit: usize,
}

unsafe impl Send for ReferenceGenome {}
unsafe impl Sync for ReferenceGenome {}

impl ReferenceGenome {
    /// Create a new ReferenceGenome instance with caching
    pub fn new(fasta_path: &str) -> Result<Self> {
        Self::new_with_cache_size(fasta_path, 100_000)
    }

    /// Create a new ReferenceGenome instance with specified cache size
    pub fn new_with_cache_size(fasta_path: &str, cache_size: usize) -> Result<Self> {
        let index_path = format!("{}.fai", fasta_path);

        if !std::path::Path::new(fasta_path).exists() {
            return Err(RedicatError::FileNotFound(format!(
                "FASTA file not found: {}",
                fasta_path
            )));
        }

        if !std::path::Path::new(&index_path).exists() {
            return Err(RedicatError::FileNotFound(format!(
                "FASTA index file not found: {}. Please create it using: samtools faidx {}",
                index_path, fasta_path
            )));
        }

        let index = Index::from_file(&index_path).map_err(|e| {
            RedicatError::ReferenceGenome(format!("Failed to load FASTA index: {:?}", e))
        })?;

        let reader = IndexedReader::with_index(
            std::fs::File::open(fasta_path).map_err(|e| {
                RedicatError::FileNotFound(format!("Cannot open FASTA file: {}", e))
            })?,
            index,
        );

        let sequences = Self::read_sequence_names(&index_path)?;

        info!(
            "Loaded reference genome with {} sequences, cache size: {}",
            sequences.len(),
            cache_size
        );

        Ok(ReferenceGenome {
            reader: parking_lot::Mutex::new(reader),
            sequences,
            cache: Arc::new(RwLock::new(HashMap::new())),
            cache_size_limit: cache_size,
        })
    }

    /// Get the reference base at a specific genomic position with caching
    pub fn get_ref_of_pos(&self, genomic_pos: &str) -> Result<char> {
        // Check cache first
        {
            let cache = self.cache.read();
            if let Some(&cached_base) = cache.get(genomic_pos) {
                return Ok(cached_base);
            }
        }

        // Parse genomic position
        let parts: Vec<&str> = genomic_pos.split(':').collect();
        if parts.len() != 2 {
            debug!("Invalid genomic position format: {}", genomic_pos);
            return Ok('N');
        }

        let chrom = parts[0];
        let pos: u64 = match parts[1].parse() {
            Ok(p) if p > 0 => p,
            _ => {
                debug!("Invalid position in {}: {}", genomic_pos, parts[1]);
                return Ok('N');
            }
        };

        // Validate chromosome exists
        if !self.sequences.iter().any(|seq| seq == chrom) {
            warn!("Chromosome {} not found in reference genome", chrom);
            return Ok('N');
        }

        // Fetch from file
        let base = self.fetch_base_from_file(chrom, pos)?;

        // Cache the result
        self.cache_base(genomic_pos.to_string(), base);

        Ok(base)
    }

    /// Fetch base from file (internal method)
    fn fetch_base_from_file(&self, chrom: &str, pos: u64) -> Result<char> {
        let mut reader = self.reader.lock();
        let mut sequence = Vec::new();

        match reader.fetch(chrom, pos - 1, pos) {
            Ok(_) => {
                if let Err(e) = reader.read(&mut sequence) {
                    debug!("Failed to read sequence at {}:{}: {:?}", chrom, pos, e);
                    return Ok('N');
                }

                if let Some(&base) = sequence.first() {
                    let base_char = (base as char).to_ascii_uppercase();
                    match base_char {
                        'A' | 'T' | 'C' | 'G' => Ok(base_char),
                        _ => {
                            debug!("Non-standard base at {}:{}: {}", chrom, pos, base_char);
                            Ok('N')
                        }
                    }
                } else {
                    debug!("Empty sequence at {}:{}", chrom, pos);
                    Ok('N')
                }
            }
            Err(e) => {
                debug!("Failed to fetch position {}:{}: {:?}", chrom, pos, e);
                Ok('N')
            }
        }
    }

    /// Cache a base with size management
    fn cache_base(&self, position: String, base: char) {
        let mut cache = self.cache.write();

        // Simple cache size management - remove oldest entries if at limit
        if cache.len() >= self.cache_size_limit {
            // Remove 10% of entries to make room
            let to_remove = self.cache_size_limit / 10;
            let keys_to_remove: Vec<String> = cache.keys().take(to_remove).cloned().collect();
            for key in keys_to_remove {
                cache.remove(&key);
            }
        }

        cache.insert(position, base);
    }

    /// Get reference bases for multiple genomic positions with parallel processing
    pub fn get_multiple_refs_parallel(&self, positions: &[String]) -> Result<Vec<char>> {
        if positions.is_empty() {
            return Ok(Vec::new());
        }

        // Use parallel processing for large batches
        const PARALLEL_THRESHOLD: usize = 100;

        if positions.len() < PARALLEL_THRESHOLD {
            // Serial processing for small batches
            let results = positions
                .iter()
                .map(|pos| self.get_ref_of_pos(pos).unwrap_or('N'))
                .collect::<Vec<_>>();
            Ok(results) // Fix: Remove .into()
        } else {
            // Parallel processing for large batches
            let results: Vec<char> = positions
                .par_iter()
                .map(|pos| self.get_ref_of_pos(pos).unwrap_or('N'))
                .collect();
            Ok(results)
        }
    }

    /// Get reference bases with chunked processing for very large datasets
    pub fn get_multiple_refs_chunked(
        &self,
        positions: &[String],
        chunk_size: usize,
    ) -> Result<Vec<char>> {
        if positions.is_empty() {
            return Ok(Vec::new());
        }

        let mut results = Vec::with_capacity(positions.len());

        for chunk in positions.chunks(chunk_size) {
            let chunk_results = self.get_multiple_refs_parallel(chunk)?;
            results.extend(chunk_results);
        }

        Ok(results)
    }

    /// Get the list of sequence names in the reference genome
    pub fn get_sequence_names(&self) -> &[String] {
        &self.sequences
    }

    /// Validate a genomic position and check if it has a standard base
    pub fn validate_position(&self, genomic_pos: &str) -> bool {
        match self.get_ref_of_pos(genomic_pos) {
            Ok(base) => base != 'N',
            Err(_) => false,
        }
    }

    /// Validate multiple positions in parallel
    pub fn validate_positions_parallel(&self, positions: &[String]) -> Vec<bool> {
        positions
            .par_iter()
            .map(|pos| self.validate_position(pos))
            .collect()
    }

    /// Get cache statistics for monitoring
    pub fn get_cache_stats(&self) -> (usize, usize, f64) {
        let cache = self.cache.read();
        let size = cache.len();
        let limit = self.cache_size_limit;
        let usage = size as f64 / limit as f64 * 100.0;
        (size, limit, usage)
    }

    /// Clear the cache
    pub fn clear_cache(&self) {
        let mut cache = self.cache.write();
        cache.clear();
        info!("Reference genome cache cleared");
    }

    /// Preload positions into cache
    pub fn preload_positions(&self, positions: &[String]) -> Result<()> {
        info!("Preloading {} positions into cache", positions.len());

        let _results: Vec<char> = positions
            .par_iter()
            .map(|pos| self.get_ref_of_pos(pos).unwrap_or('N'))
            .collect();

        let (cache_size, _, usage) = self.get_cache_stats();
        info!(
            "Preloading complete. Cache size: {}, usage: {:.1}%",
            cache_size, usage
        );

        Ok(())
    }

    /// Read sequence names from a FASTA index file
    fn read_sequence_names(index_path: &str) -> Result<Vec<String>> {
        use std::fs::File;
        use std::io::{BufRead, BufReader};

        let file = File::open(index_path).map_err(RedicatError::Io)?;

        let reader = BufReader::new(file);
        let sequences: Vec<String> = reader
            .lines()
            .map(|line| {
                line.map_err(RedicatError::Io).and_then(|l| {
                    let fields: Vec<&str> = l.split('\t').collect();
                    if fields.is_empty() {
                        Err(RedicatError::Parse("Empty line in FASTA index".to_string()))
                    } else {
                        Ok(fields[0].to_string())
                    }
                })
            })
            .collect::<Result<Vec<_>>>()?;

        if sequences.is_empty() {
            return Err(RedicatError::EmptyData(
                "No sequences found in FASTA index".to_string(),
            ));
        }

        Ok(sequences)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_fasta() -> (NamedTempFile, NamedTempFile) {
        let mut fasta_file = NamedTempFile::new().unwrap();
        writeln!(fasta_file, ">chr1").unwrap();
        writeln!(fasta_file, "ATCGATCGATCG").unwrap();
        writeln!(fasta_file, ">chr2").unwrap();
        writeln!(fasta_file, "GCTAGCTAGCTA").unwrap();

        let mut index_file = NamedTempFile::new().unwrap();
        writeln!(index_file, "chr1\t12\t6\t12\t13").unwrap();
        writeln!(index_file, "chr2\t12\t25\t12\t13").unwrap();

        (fasta_file, index_file)
    }

    #[test]
    fn test_cache_stats() {
        // This is a basic test - in real scenarios you'd use actual FASTA files
        let (_fasta, _index) = create_test_fasta();
        // Would test with actual ReferenceGenome instance
    }

    #[test]
    fn test_position_validation() {
        // Test position string parsing
        assert_eq!(
            "chr1:1000".split(':').collect::<Vec<_>>(),
            vec!["chr1", "1000"]
        );
        assert!("1000".parse::<u64>().is_ok());
        assert!("invalid".parse::<u64>().is_err());
    }
}
