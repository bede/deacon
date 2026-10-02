//! Download prebuilt indexes.

use crate::index_format::INDEX_FORMAT_VERSION;
use anyhow::{Context, Result};
use indicatif::ProgressBar;
use std::path::{Path, PathBuf};
use tracing::info;

/// Download an index to `output` or `./{name}.k{k}w{w}.idx` and return its path.
pub fn fetch(
    index_name: &str,
    kmer_length: u8,
    window_size: u8,
    output: Option<&Path>,
) -> Result<PathBuf> {
    const DEFAULT_REPOSITORY_URL: &str =
        "https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o";

    let base_url = std::env::var("DEACON_REPOSITORY_URL")
        .unwrap_or_else(|_| DEFAULT_REPOSITORY_URL.to_string());

    let filename = format!("{}.k{}w{}.idx", index_name, kmer_length, window_size);
    let url = format!("{}/deacon/{}/{}", base_url, INDEX_FORMAT_VERSION, filename);

    info!("Fetching {}", url);

    let mut response = minreq::get(&url)
        .send_lazy()
        .context("Failed to download index")?;

    if response.status_code != 200 {
        anyhow::bail!("Failed to fetch index: HTTP {}", response.status_code);
    }

    let content_length = response
        .headers
        .iter()
        .find(|(name, _)| name == "content-length")
        .and_then(|(_, value)| value.parse::<u64>().ok())
        .unwrap_or(0);

    let pb = ProgressBar::new(content_length);

    let output_path = output
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from(&filename));

    let mut temp_path = output_path.clone();
    temp_path.as_mut_os_string().push(".tmp");

    let mut file = std::fs::File::create(&temp_path).context("Failed to create temporary file")?;
    std::io::copy(&mut pb.wrap_read(&mut response), &mut file)
        .context("Failed to write index to file")?;

    std::fs::rename(&temp_path, &output_path).context("Failed to finalise index")?;

    pb.finish_and_clear();
    info!("Index saved to {}", output_path.display());

    Ok(output_path)
}
