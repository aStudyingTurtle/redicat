use anyhow::{Context, Result};
use crossbeam::channel::{bounded, Receiver};
use log::*;
use num_cpus;
use rayon::prelude::*;
use rust_htslib::bam::{IndexedReader, Read};
use rust_lapper::Lapper;
use std::{convert::TryInto, path::PathBuf, thread};

use super::intervals;
use super::types::{RegionProcessor, BYTES_IN_A_GIGABYTE, CHANNEL_SIZE_MODIFIER, CHUNKSIZE};

/// Parallel BAM/CRAM region executor driven by [`RegionProcessor`] implementations.
#[derive(Debug)]
pub struct ParGranges<R: 'static + RegionProcessor + Send + Sync> {
    reads: PathBuf,
    ref_fasta: Option<PathBuf>,
    regions_bed: Option<PathBuf>,
    regions_bcf: Option<PathBuf>,
    merge_regions: bool,
    threads: usize,
    chunksize: u32,
    channel_size_modifier: f64,
    pool: rayon::ThreadPool,
    processor: R,
}

impl<R: RegionProcessor + Send + Sync> ParGranges<R> {
    /// Create a new [`ParGranges`] executor.
    #[allow(clippy::too_many_arguments)]
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
        let requested_threads = threads.unwrap_or_else(num_cpus::get);
        let threads = std::cmp::max(requested_threads, 1);
        info!("Using {} worker threads.", threads);

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("Failed to build Rayon thread pool");

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

    /// Launch parallel processing for all configured regions.
    pub fn process(self) -> Result<Receiver<R::P>> {
        let ParGranges {
            reads,
            ref_fasta,
            regions_bed,
            regions_bcf,
            merge_regions,
            threads,
            chunksize,
            channel_size_modifier,
            pool,
            processor,
        } = self;

        let item_size = std::mem::size_of::<R::P>().max(1);
        let channel_size: usize =
            ((BYTES_IN_A_GIGABYTE as f64 * channel_size_modifier).floor() as usize / item_size)
                .saturating_mul(threads);
        info!(
            "Creating channel of length {} (* {} bytes per item)",
            channel_size, item_size
        );

        let engine = Engine {
            reads,
            ref_fasta,
            regions_bed,
            regions_bcf,
            merge_regions,
            threads,
            chunksize,
            processor,
        };

        let (sender, receiver) = bounded::<R::P>(channel_size.max(1));
        thread::spawn(move || {
            pool.install(move || {
                if let Err(err) = engine.run(sender) {
                    error!("ParGranges terminated with error: {}", err);
                }
            });
        });
        Ok(receiver)
    }
}

struct Engine<R: RegionProcessor + Send + Sync> {
    reads: PathBuf,
    ref_fasta: Option<PathBuf>,
    regions_bed: Option<PathBuf>,
    regions_bcf: Option<PathBuf>,
    merge_regions: bool,
    threads: usize,
    chunksize: u32,
    processor: R,
}

impl<R: RegionProcessor + Send + Sync> Engine<R> {
    fn run(self, sender: crossbeam::channel::Sender<R::P>) -> Result<()> {
        info!("Reading from {:?}", self.reads);
        let mut reader = IndexedReader::from_path(&self.reads)
            .with_context(|| format!("Failed to open BAM/CRAM {}", self.reads.display()))?;
        if let Err(e) = reader.set_threads(self.threads) {
            error!("Failed to set thread count to {}: {}", self.threads, e);
        }
        if let Some(ref_fasta) = &self.ref_fasta {
            reader
                .set_reference(ref_fasta)
                .with_context(|| format!("Failed to set reference {}", ref_fasta.display()))?;
        }
        let header = reader.header().to_owned();
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

        let bed_intervals = match &self.regions_bed {
            Some(path) => Some(intervals::bed_to_intervals(
                &header,
                path,
                self.merge_regions,
            )?),
            None => None,
        };
        let bcf_intervals = match &self.regions_bcf {
            Some(path) => Some(intervals::bcf_to_intervals(
                &header,
                path,
                self.merge_regions,
            )?),
            None => None,
        };

        let restricted = match (bed_intervals, bcf_intervals) {
            (Some(bed), Some(bcf)) => {
                Some(intervals::merge_intervals(bed, bcf, self.merge_regions))
            }
            (Some(bed), None) => Some(bed),
            (None, Some(bcf)) => Some(bcf),
            (None, None) => None,
        };

        let intervals = match restricted {
            Some(ivs) => ivs,
            None => intervals::header_to_intervals(&header, self.chunksize)?,
        };

        let serial_step_size = self.chunksize.saturating_mul(self.threads as u32);
        info!("Processing {} contigs", intervals.len());

        intervals
            .into_par_iter()
            .enumerate()
            .for_each(|(tid_idx, contig_intervals)| {
                let tid = tid_idx as u32;
                let (tid_end, tid_name) = match target_info.get(tid as usize) {
                    Some((len, name)) => (*len, name.clone()),
                    None => {
                        error!("Missing target info for TID {}", tid);
                        return;
                    }
                };
                if tid_end == 0 {
                    error!("Invalid target length for TID {}", tid);
                    return;
                }

                for chunk_start in (0..tid_end).step_by(serial_step_size as usize) {
                    let chunk_end = std::cmp::min(chunk_start + serial_step_size, tid_end);

                    trace!(
                        "Batch processing {}:{}-{}",
                        tid_name,
                        chunk_start,
                        chunk_end
                    );

                    let region_intervals: Vec<_> =
                        Lapper::find(&contig_intervals, chunk_start, chunk_end)
                            .map(|iv| {
                                let start = std::cmp::max(iv.start, chunk_start);
                                let stop = std::cmp::min(iv.stop, chunk_end);
                                (start, stop)
                            })
                            .collect();

                    if region_intervals.is_empty() {
                        continue;
                    }

                    region_intervals.into_par_iter().for_each_with(
                        sender.clone(),
                        |snd, (start, stop)| {
                            trace!("Processing {}:{}-{} on TID {}", tid_name, start, stop, tid);
                            let results = self.processor.process_region(tid, start, stop);
                            for item in results {
                                if snd.send(item).is_err() {
                                    warn!("Channel closed; terminating region processing early");
                                    return;
                                }
                            }
                        },
                    );
                }
            });

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::bed;
    use proptest::prelude::*;
    use rust_htslib::{bam, bcf};
    use rust_lapper::{Interval, Lapper};
    use smartstring::SmartString;
    use std::collections::{HashMap, HashSet};
    use tempfile::tempdir;

    use crate::engine::position::pileup_position::PileupPosition;
    use crate::engine::position::Position;

    struct TestProcessor;

    impl RegionProcessor for TestProcessor {
        type P = PileupPosition;

        fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P> {
            (start..stop)
                .map(|pos| {
                    let chr = SmartString::from(&tid.to_string());
                    PileupPosition::new(chr, pos)
                })
                .collect()
        }
    }

    prop_compose! {
        fn arb_iv_start(max_iv: u64)(start in 0..max_iv/2) -> u64 { start }
    }
    prop_compose! {
        fn arb_iv_size(max_iv: u64)(size in 1..max_iv/2) -> u64 { size }
    }
    prop_compose! {
        fn arb_iv(max_iv: u64)(start in arb_iv_start(max_iv), size in arb_iv_size(max_iv)) -> Interval<u64, ()> {
            Interval { start, stop: start + size, val: () }
        }
    }
    fn arb_ivs(
        max_iv: u64,
        max_ivs: usize,
    ) -> impl Strategy<Value = (Vec<Interval<u64, ()>>, u64, u64)> {
        prop::collection::vec(arb_iv(max_iv), 0..max_ivs).prop_map(|vec| {
            let mut furthest_right = 0;
            let lapper = Lapper::new(vec.clone());
            let expected = lapper.cov();
            for iv in vec.iter() {
                furthest_right = furthest_right.max(iv.stop);
            }
            (vec, expected, furthest_right)
        })
    }
    fn arb_chrs(
        max_chr: usize,
        max_iv: u64,
        max_ivs: usize,
    ) -> impl Strategy<Value = Vec<(Vec<Interval<u64, ()>>, u64, u64)>> {
        prop::collection::vec(arb_ivs(max_iv, max_ivs), 0..max_chr)
    }

    proptest! {
        #[test]
        fn interval_set(
            chromosomes in arb_chrs(4, 10_000, 1_000),
            chunksize in any::<u32>(),
            cpus in 0..num_cpus::get(),
            use_bed in any::<bool>(),
            use_vcf in any::<bool>(),
        ) {
            let tempdir = tempdir().unwrap();
            let bam_path = tempdir.path().join("test.bam");
            let bed_path = tempdir.path().join("test.bed");
            let vcf_path = tempdir.path().join("test.vcf");

            let mut header = bam::header::Header::new();
            for (i, chr) in chromosomes.iter().enumerate() {
                let mut chr_rec = bam::header::HeaderRecord::new(b"SQ");
                chr_rec.push_tag(b"SN", &i.to_string());
                chr_rec.push_tag(b"LN", &chr.2.to_string());
                header.push_record(&chr_rec);
            }
            let writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            drop(writer);
            bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

            let mut bed_writer = bed::Writer::to_file(&bed_path).unwrap();
            for (i, chr) in chromosomes.iter().enumerate() {
                for iv in chr.0.iter() {
                    let mut record = bed::Record::new();
                    record.set_start(iv.start);
                    record.set_end(iv.stop);
                    record.set_chrom(&i.to_string());
                    record.set_score(&0.to_string());
                    bed_writer.write(&record).unwrap();
                }
            }
            drop(bed_writer);

            let mut vcf_truth = HashMap::new();
            let mut vcf_header = bcf::header::Header::new();
            for (i, chr) in chromosomes.iter().enumerate() {
                vcf_header.push_record(
                    format!("##contig=<ID={},length={}>", i, chr.2).as_bytes(),
                );
            }
            let mut vcf_writer = bcf::Writer::from_path(&vcf_path, &vcf_header, true, bcf::Format::Vcf).unwrap();
            let mut record = vcf_writer.empty_record();
            for (i, chr) in chromosomes.iter().enumerate() {
                record.set_rid(Some(i as u32));
                let counter = vcf_truth.entry(i).or_insert(0);
                let mut seen = HashSet::new();
                for iv in chr.0.iter() {
                    if seen.insert(iv.start) {
                        *counter += 1;
                    }
                    record.set_pos(iv.start as i64);
                    vcf_writer.write(&record).unwrap();
                }
            }
            drop(vcf_writer);

            let par_granges_runner = ParGranges::new(
                bam_path,
                None,
                if use_bed { Some(bed_path) } else { None },
                if use_vcf { Some(vcf_path) } else { None },
                true,
                Some((cpus + 1).max(1)),
                Some(chunksize.max(1)),
                Some(0.002),
                TestProcessor,
            );
            let receiver = par_granges_runner.process().unwrap();
            let mut chrom_counts = HashMap::new();
            receiver.into_iter().for_each(|p: PileupPosition| {
                *chrom_counts.entry(p.ref_seq.parse::<usize>().unwrap()).or_insert(0u64) += 1;
            });

            for (chrom, positions) in chrom_counts.iter() {
                if use_bed && !use_vcf {
                    prop_assert_eq!(chromosomes[*chrom].1, *positions);
                } else if use_bed && use_vcf {
                    prop_assert_eq!(chromosomes[*chrom].1, *positions);
                } else if use_vcf && !use_bed {
                    prop_assert_eq!(vcf_truth.get(chrom).unwrap(), positions);
                } else {
                    prop_assert_eq!(chromosomes[*chrom].2, *positions);
                }
            }
        }
    }
}
