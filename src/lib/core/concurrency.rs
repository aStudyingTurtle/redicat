use anyhow::{Error, Result};
use log::{error, warn};

/// Set the global Rayon thread pool size to the validated value.
pub fn set_rayon_global_pools_size(size: usize) -> Result<()> {
    let cpus = determine_allowed_cpus(size)?;
    rayon::ThreadPoolBuilder::new()
        .num_threads(cpus)
        .build_global()?;
    Ok(())
}

/// Validate and normalize a requested CPU count.
pub fn determine_allowed_cpus(desired: usize) -> Result<usize> {
    if desired == 0 {
        error!("Must select > 0 threads");
        Err(Error::msg("Too few threads selected. Min 4"))
    } else if desired > num_cpus::get() {
        warn!(
            "Specified more threads than are available, using {}",
            desired
        );
        Ok(desired)
    } else {
        Ok(desired)
    }
}
