//! Download prebuilt indexes.

use crate::index_format::INDEX_FORMAT_VERSION;
use anyhow::{Context, Result};
use indicatif::ProgressBar;
use std::path::{Path, PathBuf};
use tracing::info;

/// Download an index to `output` (a file or directory) or `./{name}.k{k}w{w}.idx` and return its path.
/// `progress` shows a download bar on stderr.
pub fn fetch(
    index_name: &str,
    kmer_length: u8,
    window_size: u8,
    output: Option<&Path>,
    progress: bool,
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

    let pb = if progress {
        ProgressBar::new(content_length)
    } else {
        ProgressBar::hidden()
    };

    let output_path = resolve_output_path(output, &filename);

    let mut temp_path = output_path.clone();
    temp_path.as_mut_os_string().push(".tmp");

    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("Failed to create directory {}", parent.display()))?;
    }

    // Write to a temp file and rename once complete, removing it on failure
    let res = (|| -> Result<()> {
        let mut file =
            std::fs::File::create(&temp_path).context("Failed to create temporary file")?;
        let written = std::io::copy(&mut pb.wrap_read(&mut response), &mut file)
            .context("Failed to write index to file")?;
        if content_length > 0 && written != content_length {
            anyhow::bail!("Incomplete download: received {written} of {content_length} bytes");
        }
        std::fs::rename(&temp_path, &output_path).context("Failed to finalise index")
    })();
    pb.finish_and_clear();
    if res.is_err() {
        let _ = std::fs::remove_file(&temp_path);
    }
    res?;

    info!("Index saved to {}", output_path.display());

    Ok(output_path)
}

/// `-o` may name a file or a directory (existing, or with a trailing separator).
fn resolve_output_path(output: Option<&Path>, filename: &str) -> PathBuf {
    match output {
        Some(p) if p.is_dir() || p.to_string_lossy().ends_with(std::path::is_separator) => {
            p.join(filename)
        }
        Some(p) => p.to_path_buf(),
        None => PathBuf::from(filename),
    }
}
