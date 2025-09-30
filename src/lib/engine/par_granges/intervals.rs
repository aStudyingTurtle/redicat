use anyhow::{anyhow, Context, Result};
use bio::io::bed;
use rust_htslib::{
    bam::HeaderView,
    bcf::{Read as _, Reader},
};
use rust_lapper::{Interval, Lapper};
use std::{convert::TryInto, path::PathBuf};

pub(crate) fn header_to_intervals(
    header: &HeaderView,
    chunksize: u32,
) -> Result<Vec<Lapper<u32, ()>>> {
    let mut intervals = vec![vec![]; header.target_count() as usize];
    for tid in 0..header.target_count() {
        let tid_len: u32 = header
            .target_len(tid)
            .ok_or_else(|| anyhow!("Missing target length for TID {}", tid))?
            .try_into()
            .map_err(|_| anyhow!("Target length overflow for TID {}", tid))?;
        for start in (0..tid_len).step_by(chunksize as usize) {
            let stop = std::cmp::min(start + chunksize, tid_len);
            intervals[tid as usize].push(Interval {
                start,
                stop,
                val: (),
            });
        }
    }
    Ok(intervals.into_iter().map(Lapper::new).collect())
}

pub(crate) fn bed_to_intervals(
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

pub(crate) fn bcf_to_intervals(
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

pub(crate) fn merge_intervals(
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
