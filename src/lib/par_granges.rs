//! # ParGranges
//!
//! Iterates over chunked genomic regions in parallel.
//!
//! This module provides the core parallelization functionality for REDICAT.
//! It handles breaking up genomic regions into chunks that can be processed
//! in parallel, taking care of the complexities of overlapping regions and
//! ensuring that results are returned in the correct order.
//!
//! The main struct is [`ParGranges`], which is configured with a [`RegionProcessor`]
//! implementation that defines how each genomic region should be processed.
//!
//! # Features
//!
//! - Parallel processing of genomic regions using Rayon
//! - Support for BED and BCF/VCF files to specify regions of interest
//! - Automatic merging of overlapping intervals to avoid duplicate processing
//! - Configurable chunk size for memory management
//! - Channel-based communication between threads with configurable buffer size
//! - Support for CRAM files with reference FASTA
//!
//! # Example
//!
//! ```rust
//! use redicat_lib::par_granges::{ParGranges, RegionProcessor};
//! use std::path::PathBuf;
//!
//! struct MyProcessor;
//!
//! impl RegionProcessor for MyProcessor {
//!     type P = String; // Your output type
//!
//!     fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P> {
//!         // Process the region and return results
//!         vec![format!("Processed region {}:{}-{}", tid, start, stop)]
//!     }
//! }
//!
//! let processor = MyProcessor;
//! let par_granges_runner = ParGranges::new(
//!     PathBuf::from("test.bam"),
//!     None, // ref_fasta
//!     None, // regions_bcf
//!     None, // regions_bed
//!     true, // merge_regions
//!     None, // threads
//!     None, // chunksize
//!     None, // channel_size_modifier
//!     processor,
//! );
//! ```

use anyhow::{anyhow, Context, Result};
use bio::io::bed;
use crossbeam::channel::{bounded, Receiver};
use lazy_static::lazy_static;
use log::*;
use num_cpus;
use rayon::prelude::*;
use rust_htslib::{
    bam::{HeaderView, IndexedReader, Read},
    bcf::{Read as bcfRead, Reader},
};
use rust_lapper::{Interval, Lapper};
use serde::Serialize;
use std::{convert::TryInto, path::PathBuf, thread};

const BYTES_INA_GIGABYTE: usize = 1024 * 1024 * 1024;

/// A modifier to apply to the channel size formular that is (BYTES_INA_GIGABYTE * channel_size_modifier) * threads / size_of(R::P)
/// 0.15 roughly corresponds to 1_000_000 PileupPosition objects per thread with some wiggle room.
pub const CHANNEL_SIZE_MODIFIER: f64 = 0.5;

/// The ideal number of basepairs each worker will receive. Total bp in memory at one time = `threads` * `chunksize`
pub const CHUNKSIZE: u32 = 500_000;

lazy_static! {
    /// CHANNEL_SIZE_MODIFIER as a str
    pub static ref CHANNEL_SIZE_MODIFIER_STR: String = CHANNEL_SIZE_MODIFIER.to_string();

    /// CHUNKSIZE as a str
    pub static ref CHUNKSIZE_STR: String = CHUNKSIZE.to_string();
}

/// RegionProcessor defines the methods that must be implemented to process a region
///
/// Implementors of this trait define how genomic regions should be processed.
/// The associated type `P` represents the output of processing a region.
pub trait RegionProcessor {
    /// The type of values returned by [`process_region`].
    ///
    /// This type must implement `Serialize` so that results can be sent
    /// between threads, and `Send` and `Sync` for thread safety.
    ///
    /// [`process_region`]: #method.process_region
    type P: 'static + Send + Sync + Serialize;

    /// Process a genomic region defined by tid, start, and stop.
    ///
    /// This method is called for each chunk of genomic regions that needs
    /// to be processed. The implementation should:
    ///
    /// 1. Open the BAM/CRAM file
    /// 2. Fetch the specified region
    /// 3. Process the reads/alignments in that region
    /// 4. Return the results as a Vec of type `P`
    ///
    /// # Arguments
    ///
    /// * `tid` - The target ID (chromosome) to process
    /// * `start` - The start position (0-based)
    /// * `stop` - The stop position (0-based, exclusive)
    ///
    /// # Returns
    ///
    /// A vector of processed results of type `P`
    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P>;
}

/// ParGranges holds all the information and configuration needed to launch the
/// [`ParGranges::process`].
///
/// [`ParGranges::process`]: #method.process
#[derive(Debug)]
pub struct ParGranges<R: 'static + RegionProcessor + Send + Sync> {
    /// Path to an indexed BAM / CRAM file
    reads: PathBuf,
    /// Optional reference file for CRAM
    ref_fasta: Option<PathBuf>,
    /// Optional path to a BED file to restrict the regions iterated over
    regions_bed: Option<PathBuf>,
    /// Optional path to a BCF/VCF file to restrict the regions iterated over
    regions_bcf: Option<PathBuf>,
    /// If `regions_bed` and or `regions_bcf` is specified, and this is true, merge any overlapping regions to avoid duplicate output.
    merge_regions: bool,
    /// Number of threads this is allowed to use, uses all if None
    threads: usize,
    /// The ideal number of basepairs each worker will receive. Total bp in memory at one time = `threads` * `chunksize`
    chunksize: u32,
    /// A modifier to apply to the channel size formular that is (BYTES_INA_GIGABYTE * channel_size_modifier) * threads / size_of(R::P)
    channel_size_modifier: f64,
    /// The rayon threadpool to operate in
    pool: rayon::ThreadPool,
    /// The implementation of [RegionProcessor] that will be used to process regions
    processor: R,
}

impl<R: RegionProcessor + Send + Sync> ParGranges<R> {
    /// Create a ParIO object
    ///
    /// # Arguments
    ///
    /// * `reads`- path to an indexed BAM/CRAM
    /// * `ref_fasta`- path to an indexed reference file for CRAM
    /// * `regions_bed`- Optional BED file path restricting the regions to be examined
    /// * `regions_bcf`- Optional BCF/VCF file path restricting the regions to be examined
    /// * `merge_regions` - If `regions_bed` and or `regions_bcf` is specified, and this is true, merge any overlapping regions to avoid duplicate output.
    /// * `threads`- Optional threads to restrict the number of threads this process will use, defaults to all
    /// * `chunksize`- optional argument to change the default chunksize of 1_000_000. `chunksize` determines the number of bases
    ///                each worker will get to work on at one time.
    /// * `channel_size_modifier`- Optional argument to modify the default size ration of the channel that `R::P` is sent on.
    ///                 formula is: ((BYTES_INA_GIGABYTE * channel_size_modifier) * threads) / size_of(R::P)
    /// * `processor`- Something that implements [`RegionProcessor`]
    pub fn new(
        reads: PathBuf,
        ref_fasta: Option<PathBuf>,
        regions_bed: Option<PathBuf>,
        regions_bcf: Option<PathBuf>,
        merge_regions: bool,
        threads: Option<usize>,
        chunksize: Option<u32>,
        channel_size_modifier: Option<f64>,
        processor: R,
    ) -> Self {
        let threads = if let Some(threads) = threads {
            threads
        } else {
            num_cpus::get()
        };

        // Keep two around for main thread and thread running the pool
        let threads = std::cmp::max(threads, 1);
        info!("Using {} threads.", threads);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();

        info!("Using {} worker threads.", threads);
        Self {
            reads,
            ref_fasta,
            regions_bed,
            regions_bcf,
            merge_regions,
            threads,
            chunksize: chunksize.unwrap_or(CHUNKSIZE),
            channel_size_modifier: channel_size_modifier.unwrap_or(CHANNEL_SIZE_MODIFIER),
            pool,
            processor,
        }
    }

    /// Process each region.
    ///
    /// This method splits the sequences in the BAM/CRAM header into `chunksize` * `self.threads` regions (aka 'super chunks').
    /// It then queries that 'super chunk' against the intervals (either the BED file, or the whole genome broken up into `chunksize`
    /// regions). The results of that query are then processed by a pool of workers that apply `process_region` to reach interval to
    /// do  analysis on. The collected result for each region is then sent back over the returned `Receiver<R::P>` channel
    /// for the caller to use. The results will be returned in order according to the order of the intervals used to drive this method.
    ///
    /// While one 'super chunk' is being worked on by all workers, the last 'super chunks' results are being printed to either to
    /// a file or to STDOUT, in order.
    ///
    /// Note, a common use case of this will be to fetch a region and do a pileup. The bounds of bases being looked at should still be
    /// checked since a fetch will pull all reads that overlap the region in question.
    pub fn process(self) -> Result<Receiver<R::P>> {
        let channel_size: usize = ((BYTES_INA_GIGABYTE as f64 * self.channel_size_modifier).floor()
            as usize
            / std::mem::size_of::<R::P>())
        .saturating_mul(self.threads);
        info!(
            "Creating channel of length {:?} (* {} bytes to get mem)",
            channel_size,
            std::mem::size_of::<R::P>()
        );
        let (snd, rxv) = bounded::<R::P>(channel_size);
        thread::spawn(move || {
            self.pool.install(|| {
                info!("Reading from {:?}", self.reads);
                let mut reader = match IndexedReader::from_path(&self.reads) {
                    Ok(r) => r,
                    Err(e) => {
                        log::error!(
                            "Failed to open BAM/CRAM file '{}': {}",
                            self.reads.display(),
                            e
                        );
                        return;
                    }
                };
                // reader
                //     .set_threads(self.threads)
                //     .expect("Failed to set number of threads");
                if let Err(e) = reader.set_threads(self.threads) {
                    log::error!("Failed to set thread count to {}: {}", self.threads, e);
                    return;
                }
                // If passed add ref_fasta
                if let Some(ref_fasta) = &self.ref_fasta {
                    if let Err(e) = reader.set_reference(ref_fasta) {
                        log::error!(
                            "Failed to set reference file '{}': {}",
                            ref_fasta.display(),
                            e
                        );
                        return;
                    }
                }
                // Get a copy of the header and extract necessary data upfront
                let header = reader.header().to_owned();

                // Pre-extract target lengths and names to avoid sharing HeaderView across threads
                let target_info: Vec<(u32, String)> = (0..header.target_count())
                    .map(|tid| {
                        let len = header
                            .target_len(tid)
                            .and_then(|len| len.try_into().ok())
                            .unwrap_or(0);
                        let name = std::str::from_utf8(header.tid2name(tid))
                            .unwrap_or("unknown")
                            .to_string();
                        (len, name)
                    })
                    .collect();

                // Work out if we are restricted to a subset of sites
                let bed_intervals = if let Some(regions_bed) = &self.regions_bed {
                    Some(
                        Self::bed_to_intervals(&header, regions_bed, self.merge_regions)
                            .expect("Parsed BED to intervals"),
                    )
                } else {
                    None
                };
                let bcf_intervals = if let Some(regions_bcf) = &self.regions_bcf {
                    Some(
                        Self::bcf_to_intervals(&header, regions_bcf, self.merge_regions)
                            .expect("Parsed BCF/VCF to intervals"),
                    )
                } else {
                    None
                };
                let restricted_ivs = match (bed_intervals, bcf_intervals) {
                    (Some(bed_ivs), Some(bcf_ivs)) => {
                        Some(Self::merge_intervals(bed_ivs, bcf_ivs, self.merge_regions))
                    }
                    (Some(bed_ivs), None) => Some(bed_ivs),
                    (None, Some(bcf_ivs)) => Some(bcf_ivs),
                    (None, None) => None,
                };

                let intervals = if let Some(restricted) = restricted_ivs {
                    restricted
                } else {
                    Self::header_to_intervals(&header, self.chunksize)
                        .expect("Parsed BAM/CRAM header to intervals")
                };

                // The number positions to try to process in one batch
                let serial_step_size = self
                    .chunksize
                    .checked_mul(self.threads as u32)
                    .unwrap_or(u32::MAX); // aka superchunk
                info!("Processing {} intervals", intervals.len());
                intervals
                    .into_par_iter()
                    .enumerate()
                    .for_each(|(tid, intervals)| {
                        let tid: u32 = tid as u32;
                        let (tid_end, tid_name) = match target_info.get(tid as usize) {
                            Some((len, name)) => (*len, name.clone()),
                            None => {
                                log::error!("Missing target info for TID {}", tid);
                                return;
                            }
                        };
                        if tid_end == 0 {
                            log::error!("Invalid target length for TID {}", tid);
                            return;
                        }
                        info!("Processing {}:0-{}", tid_name, tid_end);

                        for chunk_start in (0..tid_end).step_by(serial_step_size as usize) {
                            let chunk_end =
                                std::cmp::min(chunk_start as u32 + serial_step_size, tid_end);

                            trace!(
                                "Batch Processing {}:{}-{}",
                                tid_name,
                                chunk_start,
                                chunk_end
                            );

                            // Pre-collect intervals for parallel processing with memory trade-off
                            let region_intervals: Vec<_> =
                                Lapper::<u32, ()>::find(&intervals, chunk_start, chunk_end)
                                    .map(|iv| {
                                        let start = std::cmp::max(iv.start, chunk_start);
                                        let stop = std::cmp::min(iv.stop, chunk_end);
                                        (start, stop)
                                    })
                                    .collect();

                            if !region_intervals.is_empty() {
                                // Use parallel processing with thread-local buffers for maximum throughput
                                region_intervals.into_par_iter().for_each_with(
                                    snd.clone(),
                                    |snd, (start, stop)| {
                                        trace!(
                                            "Processing {}:{}-{} via optimized parallel",
                                            tid_name,
                                            start,
                                            stop
                                        );
                                        let results =
                                            self.processor.process_region(tid, start, stop);

                                        // Processing results
                                        for item in results {
                                            if snd.send(item).is_err() {
                                                log::warn!(
                                                    "Channel closed, writer may have disconnected"
                                                );
                                                return;
                                            }
                                        }
                                    },
                                );
                            }
                        }
                    });
            });
        });
        Ok(rxv)
    }

    // Convert the header into intervals of equally sized chunks. The last interval may be short.
    fn header_to_intervals(header: &HeaderView, chunksize: u32) -> Result<Vec<Lapper<u32, ()>>> {
        let mut intervals = vec![vec![]; header.target_count() as usize];
        for tid in 0..(header.target_count()) {
            let tid_len: u32 = header
                .target_len(tid)
                .ok_or_else(|| anyhow!("Missing target length for TID {}", tid))?
                .try_into()
                .map_err(|_| anyhow!("Target length overflow for TID {}", tid))?;
            for start in (0..tid_len).step_by(chunksize as usize) {
                let stop = std::cmp::min(start as u32 + chunksize, tid_len);
                intervals[tid as usize].push(Interval {
                    start: start as u32,
                    stop: stop,
                    val: (),
                });
            }
        }
        Ok(intervals.into_iter().map(|ivs| Lapper::new(ivs)).collect())
    }

    /// Read a bed file into a vector of lappers with the index representing the TID
    /// if `merge' is true then any overlapping intervals in the sets will be merged.
    // TODO add a proper error message
    fn bed_to_intervals(
        header: &HeaderView,
        bed_file: &PathBuf,
        merge: bool,
    ) -> Result<Vec<Lapper<u32, ()>>> {
        let mut bed_reader = bed::Reader::from_file(bed_file)?;
        let mut intervals = vec![vec![]; header.target_count() as usize];
        for (i, record) in bed_reader.records().enumerate() {
            let record = record?;
            let tid = header.tid(record.chrom().as_bytes()).ok_or_else(|| {
                anyhow!(
                    "Chromosome '{}' not found in BAM/CRAM header",
                    record.chrom()
                )
            })?;
            let start = record
                .start()
                .try_into()
                .with_context(|| format!("BED record {} is invalid: unable to parse start", i))?;
            let stop = record
                .end()
                .try_into()
                .with_context(|| format!("BED record {} is invalid: unable to parse stop", i))?;
            if stop < start {
                return Err(anyhow!("BED record {} is invalid: stop < start", i));
            }
            intervals[tid as usize].push(Interval {
                start,
                stop,
                val: (),
            });
        }

        Ok(intervals
            .into_iter()
            .map(|ivs| {
                let mut lapper = Lapper::new(ivs);
                if merge {
                    lapper.merge_overlaps();
                }
                lapper
            })
            .collect())
    }

    /// Read a BCF/VCF file into a vector of lappers with index representing the TID
    /// if `merge' is true then any overlapping intervals in the sets will be merged.
    fn bcf_to_intervals(
        header: &HeaderView,
        bcf_file: &PathBuf,
        merge: bool,
    ) -> Result<Vec<Lapper<u32, ()>>> {
        let mut bcf_reader = Reader::from_path(bcf_file)
            .map_err(|e| anyhow!("Error opening BCF/VCF file '{}': {}", bcf_file.display(), e))?;
        let bcf_header_reader = Reader::from_path(bcf_file)
            .map_err(|e| anyhow!("Error opening BCF/VCF file '{}': {}", bcf_file.display(), e))?;
        let bcf_header = bcf_header_reader.header();
        let mut intervals = vec![vec![]; header.target_count() as usize];
        // TODO: validate the headers against eachother
        for record in bcf_reader.records() {
            let record = record?;
            let record_rid = bcf_header
                .rid2name(record.rid().unwrap())
                .map_err(|e| anyhow!("Error getting record RID name: {}", e))?;
            let tid = header.tid(record_rid).ok_or_else(|| {
                anyhow!(
                    "Chromosome '{}' not found in BAM/CRAM header",
                    String::from_utf8_lossy(record_rid)
                )
            })?;
            let pos: u32 = record
                .pos()
                .try_into()
                .map_err(|_| anyhow!("Got a negative value for pos"))?;
            intervals[tid as usize].push(Interval {
                start: pos,
                stop: pos + 1,
                val: (),
            });
        }

        Ok(intervals
            .into_iter()
            .map(|ivs| {
                let mut lapper = Lapper::new(ivs);
                if merge {
                    lapper.merge_overlaps();
                }
                lapper
            })
            .collect())
    }

    /// Merge two sets of restriction intervals together
    /// if `merge' is true then any overlapping intervals in the sets will be merged.
    fn merge_intervals(
        a_ivs: Vec<Lapper<u32, ()>>,
        b_ivs: Vec<Lapper<u32, ()>>,
        merge: bool,
    ) -> Vec<Lapper<u32, ()>> {
        let mut intervals = vec![vec![]; a_ivs.len()];
        for (i, (a_lapper, b_lapper)) in a_ivs.into_iter().zip(b_ivs.into_iter()).enumerate() {
            intervals[i] = a_lapper.into_iter().chain(b_lapper.into_iter()).collect();
        }
        intervals
            .into_iter()
            .map(|ivs| {
                let mut lapper = Lapper::new(ivs);
                if merge {
                    lapper.merge_overlaps();
                }
                lapper
            })
            .collect()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use bio::io::bed;
    use num_cpus;
    use proptest::prelude::*;
    use rust_htslib::{bam, bcf};
    use rust_lapper::{Interval, Lapper};
    use std::collections::{HashMap, HashSet};
    use tempfile::tempdir;
    // The purpose of these tests is to demonstrate that positions are covered once under a variety of circumstances

    prop_compose! {
        fn arb_iv_start(max_iv: u64)(start in 0..max_iv/2) -> u64 { start }
    }
    prop_compose! {
        fn arb_iv_size(max_iv: u64)(size in 1..max_iv/2) -> u64 { size }
    }
    prop_compose! {
        // Create an arbitrary interval where the min size == max_iv / 2
        fn arb_iv(max_iv: u64)(start in arb_iv_start(max_iv), size in arb_iv_size(max_iv)) -> Interval<u64, ()> {
            Interval {start, stop: start + size, val: ()}
        }
    }
    // Create an arbitrary number of intervals along with the expected number of positions they cover
    fn arb_ivs(
        max_iv: u64,    // max iv size
        max_ivs: usize, // max number of intervals
    ) -> impl Strategy<Value = (Vec<Interval<u64, ()>>, u64, u64)> {
        prop::collection::vec(arb_iv(max_iv), 0..max_ivs).prop_map(|vec| {
            let mut furthest_right = 0;
            let lapper = Lapper::new(vec.clone());
            let expected = lapper.cov();
            for iv in vec.iter() {
                if iv.stop > furthest_right {
                    furthest_right = iv.stop;
                }
            }
            (vec, expected, furthest_right)
        })
    }
    // Create arbitrary number of contigs with arbitrary intervals each
    fn arb_chrs(
        max_chr: usize, // number of chromosomes to use
        max_iv: u64,    // max interval size
        max_ivs: usize, // max number of intervals
    ) -> impl Strategy<Value = Vec<(Vec<Interval<u64, ()>>, u64, u64)>> {
        prop::collection::vec(arb_ivs(max_iv, max_ivs), 0..max_chr)
    }
    // An empty BAM with correct header
    // A BED file with the randomly generated intervals (with expected number of positions)
    // proptest generate random chunksize, cpus
    proptest! {
        #[test]
        // add random chunksize and random cpus
        // NB: using any larger numbers for this tends to blow up the test runtime
        fn interval_set(chromosomes in arb_chrs(4, 10_000, 1_000), chunksize in any::<u32>(), cpus in 0..num_cpus::get(), use_bed in any::<bool>(), use_vcf in any::<bool>()) {
            let tempdir = tempdir().unwrap();
            let bam_path = tempdir.path().join("test.bam");
            let bed_path = tempdir.path().join("test.bed");
            let vcf_path = tempdir.path().join("test.vcf");

            // Build a BAM
            let mut header = bam::header::Header::new();
            for (i,chr) in chromosomes.iter().enumerate() {
                let mut chr_rec = bam::header::HeaderRecord::new(b"SQ");
                chr_rec.push_tag(b"SN", &i.to_string());
                chr_rec.push_tag(b"LN", &chr.2.to_string()); // set len as max observed
                header.push_record(&chr_rec);
            }
            let writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).expect("Opened test.bam for writing");
            drop(writer); // force flush the writer so the header info is written
            bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

            // Build a bed
            let mut writer = bed::Writer::to_file(&bed_path).expect("Opened test.bed for writing");
            for (i, chr) in chromosomes.iter().enumerate() {
                for iv in chr.0.iter() {
                    let mut record = bed::Record::new();
                    record.set_start(iv.start);
                    record.set_end(iv.stop);
                    record.set_chrom(&i.to_string());
                    record.set_score(&0.to_string());
                    writer.write(&record).expect("Wrote to test.bed");
                }
            }
            drop(writer); // force flush

            // Build a VCF file
            let mut vcf_truth = HashMap::new();
            let mut header = bcf::header::Header::new();
            for (i,chr) in chromosomes.iter().enumerate() {
                header.push_record(format!("##contig=<ID={},length={}>", &i.to_string(), &chr.2.to_string()).as_bytes());
            }
            let mut writer = bcf::Writer::from_path(&vcf_path, &header, true, bcf::Format::Vcf).expect("Failed to open test.vcf for writing");
            let mut record = writer.empty_record();
            for (i, chr) in chromosomes.iter().enumerate() {
                record.set_rid(Some(i as u32));
                let counter = vcf_truth.entry(i).or_insert(0);
                let mut seen = HashSet::new();
                for iv in chr.0.iter() {
                    if !seen.contains(&iv.start) {
                        *counter += 1;
                        seen.insert(iv.start);
                    }
                    record.set_pos(iv.start as i64);
                    writer.write(&record).expect("Failed to write to test.vcf")
                }
            }

            drop(writer); // force flush
            // Create the processor with a dumb impl of processing that just returns positions with no counting
            let test_processor = TestProcessor {};
            let par_granges_runner = ParGranges::new(
                bam_path,
                None,
                if use_bed { Some(bed_path) } else { None }, // do one with regions
                if use_vcf { Some(vcf_path) } else { None }, // do one with vcf regions
                true,
                Some(cpus),
                Some(chunksize),
                Some(0.002),
                test_processor
            );
            let receiver = par_granges_runner.process().expect("Launch ParGranges Process");
            let mut chrom_counts = HashMap::new();
            receiver.into_iter().for_each(|p: PileupPosition| {
                let positions = chrom_counts.entry(p.ref_seq.parse::<usize>().expect("parsed chr")).or_insert(0u64);
                *positions += 1
            });

            // Validate that for each chr we get the expected number of bases
            for (chrom, positions) in chrom_counts.iter() {
                if use_bed  && !use_vcf {
                    // if this was with bed, should be equal to .1
                    prop_assert_eq!(chromosomes[*chrom].1, *positions, "chr: {}, expected: {}, found: {}", chrom, chromosomes[*chrom].1, positions);
                } else if use_bed && use_vcf {
                    // if this was with bed, should be equal to .1, bed restrictions and vcf restrctions should overlap
                    prop_assert_eq!(chromosomes[*chrom].1, *positions, "chr: {}, expected: {}, found: {}", chrom, chromosomes[*chrom].1, positions);
                } else if use_vcf && !use_bed {
                    // total positions should be equal to the number of records for that chr in the vcf
                    prop_assert_eq!(vcf_truth.get(chrom).unwrap(), positions, "chr: {}, expected: {}, found: {}", chrom, chromosomes[*chrom].1, positions);
                } else {
                    // if this was bam only, should be equal to rightmost postion
                    prop_assert_eq!(chromosomes[*chrom].2, *positions, "chr: {}, expected: {}, found: {}", chrom, chromosomes[*chrom].2, positions);
                }
            }

        }
    }

    use crate::position::{pileup_position::PileupPosition, Position};
    use smartstring::SmartString;
    use std::path::Path;
    struct TestProcessor {}
    impl RegionProcessor for TestProcessor {
        type P = PileupPosition;

        fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P> {
            let mut results = vec![];
            for i in start..stop {
                let chr = SmartString::from(&tid.to_string());
                let pos = PileupPosition::new(chr, i);
                results.push(pos);
            }
            results
        }
    }

    #[test]
    fn test_chr22_bam() {
        let bam_path = Path::new("test/chr22.bam");
        // Only run this test if the file exists
        if bam_path.exists() {
            // Create a simple processor that just counts positions
            let test_processor = TestProcessor {};

            // Create a ParGranges runner with a small chunksize for testing
            let par_granges_runner = ParGranges::new(
                bam_path.to_path_buf(),
                None,          // ref_fasta
                None,          // regions_bed
                None,          // regions_bcf
                true,          // merge_regions
                Some(2),       // threads
                Some(100_000), // chunksize - smaller for testing
                Some(0.5),     // channel_size_modifier
                test_processor,
            );

            // Process and collect some results
            let receiver = par_granges_runner
                .process()
                .expect("Launch ParGranges Process");

            // Collect first 100 results to verify it works (smaller number to avoid channel issues)
            let mut count = 0;
            for _ in receiver.into_iter().take(100) {
                count += 1;
            }
            info!("Collected {} results", count);
            // Verify we got some results
            assert_eq!(count, 100, "Should have processed 100 positions");
        }
        // If file doesn't exist, test silently passes
    }
}
