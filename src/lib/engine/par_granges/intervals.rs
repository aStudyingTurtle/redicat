use anyhow::{anyhow, Context, Result};
use bio::io::bed;
use noodles::{bcf, sam, vcf};
use rust_lapper::{Interval, Lapper};
use std::{convert::{TryFrom, TryInto}, fs::File, path::PathBuf};

pub(crate) fn header_to_intervals(
    header: &sam::Header,
    chunksize: u32,
) -> Result<Vec<Lapper<u32, ()>>> {
    let reference_sequences = header.reference_sequences();
    let mut intervals = vec![Vec::new(); reference_sequences.len()];

    for (tid, (_name, sequence)) in reference_sequences.iter().enumerate() {
        let length: u32 = sequence
            .length()
            .get()
            .try_into()
            .map_err(|_| anyhow!("Target length overflow for TID {}", tid))?;

        for start in (0..length).step_by(chunksize as usize) {
            let stop = std::cmp::min(start + chunksize, length);
            intervals[tid].push(Interval {
                start,
                stop,
                val: (),
            });
        }
    }
    Ok(intervals.into_iter().map(Lapper::new).collect())
}

pub(crate) fn bed_to_intervals(
    header: &sam::Header,
    bed_file: &PathBuf,
    merge: bool,
) -> Result<Vec<Lapper<u32, ()>>> {
    let mut bed_reader = bed::Reader::from_file(bed_file)?;
    let reference_sequences = header.reference_sequences();
    let mut intervals = vec![Vec::new(); reference_sequences.len()];
    for (i, record) in bed_reader.records().enumerate() {
        let record = record?;
        let tid = reference_sequences
            .get_index_of(record.chrom().as_bytes())
            .ok_or_else(|| {
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
    header: &sam::Header,
    bcf_file: &PathBuf,
    merge: bool,
) -> Result<Vec<Lapper<u32, ()>>> {
    let file = File::open(bcf_file)
        .with_context(|| format!("Error opening BCF/VCF file '{}'", bcf_file.display()))?;
    let mut bcf_reader = bcf::io::Reader::new(file);
    let vcf_header = bcf_reader
        .read_header()
        .with_context(|| format!("Error reading BCF/VCF header from '{}'", bcf_file.display()))?;
    let string_maps = vcf::header::StringMaps::try_from(&vcf_header)
        .map_err(|e| anyhow!("Error building string maps: {}", e))?;

    let reference_sequences = header.reference_sequences();
    let mut intervals = vec![Vec::new(); reference_sequences.len()];

    for record_result in bcf_reader.records() {
        let record = record_result
            .with_context(|| format!("Error reading record from '{}'", bcf_file.display()))?;

        let contig = record
            .reference_sequence_name(&string_maps)
            .map_err(|e| anyhow!("Error getting record reference sequence name: {}", e))?;

        let tid = reference_sequences
            .get_index_of(contig.as_bytes())
            .ok_or_else(|| {
                anyhow!(
                    "Chromosome '{}' not found in BAM/CRAM header",
                    contig
                )
            })?;

        let position = record
            .variant_start()
            .transpose()
            .map_err(|e| anyhow!("Error parsing variant start: {}", e))?
            .ok_or_else(|| anyhow!("BCF record is missing a variant start position"))?;

        let start = usize::from(position)
            .checked_sub(1)
            .ok_or_else(|| anyhow!("Got an invalid value for pos"))?;
        let pos: u32 = start
            .try_into()
            .map_err(|_| anyhow!("Position overflow for '{}'", contig))?;

        intervals[tid].push(Interval {
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
