mod args;
mod processor;

use anyhow::Result;
use log::info;
use redicat_lib::core::read_filter::DefaultReadFilter;
use redicat_lib::engine::{par_granges, ParGranges};
use redicat_lib::utils;
use std::sync::Arc;

use crate::commands::{common, is_standard_contig};

pub use args::{BulkArgs, BulkConfig};
use processor::BaseProcessor;

/// Execute the `bulk` command end-to-end.
pub fn run_bulk(args: BulkArgs) -> Result<()> {
    let config: BulkConfig = args.into();

    info!("Running redicat bulk on {:?}", config.reads);
    let threads = utils::determine_allowed_cpus(config.threads)?;

    let output_path = common::ensure_gz_path(&config.output);
    utils::make_parent_dirs(&output_path)?;

    let mut writer = utils::get_writer(&Some(output_path.clone()), true, true, 1, 6)?;

    let read_filter = Arc::new(DefaultReadFilter::new(config.mapquality));

    let allowed_tids = if config.all_contigs {
        None
    } else {
        Some(common::collect_tids_with_filter(&config.reads, |name| {
            is_standard_contig(name)
        })?)
    };

    let processor = BaseProcessor::from_config(&config, read_filter, allowed_tids)?;

    let runner = ParGranges::new(
        config.reads.clone(),
        None,
        None,
        None,
        false,
        Some(threads),
        Some(config.chunksize),
        Some(par_granges::CHANNEL_SIZE_MODIFIER),
        processor,
    );

    let receiver = runner.process()?;
    for pos in receiver.into_iter() {
        writer.serialize(pos)?;
    }

    writer.flush()?;
    info!("Bulk processing complete -> {:?}", output_path);
    Ok(())
}
