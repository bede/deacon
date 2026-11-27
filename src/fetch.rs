//! List and fetch prebuilt indexes

mod manifest;

use anyhow::{Context, Result, bail};
use aws_lc_rs::digest::{self, SHA256};
use indicatif::ProgressBar;
use manifest::{Asset, Manifest};
use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use tracing::info;

const DEFAULT_MANIFEST_URL: &str = "https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/3/manifest.json";
const MAX_MANIFEST_BYTES: u64 = 4 * 1024 * 1024;

/// Write the remote index listing
pub fn list(output: impl Write) -> Result<()> {
    let (manifest, _) = load_manifest()?;
    manifest.write_listing(output)
}

/// Fetch an exact index by name and k/w, reusing local files unless `force`
pub fn fetch(
    index_name: &str,
    kmer_length: u8,
    window_size: u8,
    output: Option<&Path>,
    force: bool,
    progress: bool,
) -> Result<PathBuf> {
    let filename = format!("{index_name}.k{kmer_length}w{window_size}.idx");
    manifest::validate_filename(&filename)?;
    let output_path = resolve_output_path(output, &filename);
    if !force && output_path.is_file() {
        info!(
            "Index exists at {}, skipping download (--force to replace)",
            output_path.display()
        );
        return Ok(output_path);
    }
    fetch_named(Some(&filename), output, force, progress)
}

/// Fetch by filename, group or alias (None for the catalog default).
/// Verify size and SHA-256. Output accepts a file or directory (default: asset filename).
pub fn fetch_named(
    name: Option<&str>,
    output: Option<&Path>,
    force: bool,
    progress: bool,
) -> Result<PathBuf> {
    let (manifest, directory) = load_manifest()?;
    let name = match name {
        Some(name) => name,
        None => {
            let name = manifest.default_index.as_deref().context(
                "The manifest has no default index. Run 'deacon index list' and pass a name",
            )?;
            info!("Using default index {name}");
            name
        }
    };
    let (filename, asset) = manifest.resolve(name)?;
    let output_path = resolve_output_path(output, filename);

    if !force && output_path.try_exists()? {
        verify_existing(&output_path, asset).with_context(|| {
            format!(
                "Existing file {} could not be verified (--force to replace)",
                output_path.display()
            )
        })?;
        info!("Verified {}, skipping download", output_path.display());
        return Ok(output_path);
    }

    let parent = output_path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(parent)
        .with_context(|| format!("Failed to create directory {}", parent.display()))?;
    let mut staging = tempfile::Builder::new();
    staging.prefix(".deacon-");
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        staging.permissions(fs::Permissions::from_mode(0o666));
    }
    let mut temporary = staging
        .tempfile_in(parent)
        .context("Failed to create temporary file")?;

    let url = asset_url(&directory, filename);
    info!("Fetching {filename}");
    let response = minreq::get(&url)
        .send_lazy()
        .context("Failed to download index")?;
    if response.status_code != 200 {
        bail!("Failed to fetch {filename}: HTTP {}", response.status_code);
    }

    let bar = if progress {
        ProgressBar::new(asset.size_bytes)
    } else {
        ProgressBar::hidden()
    };
    let result = copy_and_verify(bar.wrap_read(response), temporary.as_file_mut(), asset);
    bar.finish_and_clear();
    result.with_context(|| format!("Failed to verify {filename}"))?;

    if force {
        temporary.persist(&output_path)
    } else {
        temporary.persist_noclobber(&output_path)
    }
    .with_context(|| format!("Failed to save index to {}", output_path.display()))?;

    info!("Index saved to {}", output_path.display());
    Ok(output_path)
}

fn load_manifest() -> Result<(Manifest, String)> {
    let url = match std::env::var("DEACON_REPOSITORY_URL") {
        Ok(url) => url,
        Err(std::env::VarError::NotPresent) => DEFAULT_MANIFEST_URL.to_owned(),
        Err(err) => return Err(err).context("Invalid DEACON_REPOSITORY_URL"),
    };
    manifest_directory(&url)?;
    let response = minreq::get(&url)
        .with_timeout(15)
        .send_lazy()
        .context("Failed to fetch index manifest")?;
    if response.status_code != 200 {
        bail!(
            "Failed to fetch index manifest: HTTP {}",
            response.status_code
        );
    }
    let directory = manifest_directory(&response.url)?.to_owned();
    let mut body = Vec::new();
    response
        .take(MAX_MANIFEST_BYTES + 1)
        .read_to_end(&mut body)
        .context("Failed to read index manifest")?;
    if body.len() as u64 > MAX_MANIFEST_BYTES {
        bail!("Index manifest exceeds {MAX_MANIFEST_BYTES} bytes");
    }
    Ok((Manifest::parse(&body)?, directory))
}

fn manifest_directory(url: &str) -> Result<&str> {
    let url = url.split(['?', '#']).next().unwrap_or_default();
    let Some((scheme, rest)) = url.split_once("://") else {
        bail!("Expected a full HTTP(S) index manifest URL");
    };
    let Some((directory, filename)) = url.rsplit_once('/') else {
        bail!("Expected a full HTTP(S) index manifest URL");
    };
    if !matches!(scheme, "http" | "https")
        || !rest
            .split_once('/')
            .is_some_and(|(host, _)| !host.is_empty())
        || filename.is_empty()
    {
        bail!("Expected a full HTTP(S) index manifest URL");
    }
    Ok(directory)
}

fn asset_url(directory: &str, filename: &str) -> String {
    let mut url = format!("{directory}/");
    for byte in filename.bytes() {
        if byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b'~') {
            url.push(char::from(byte));
        } else {
            write!(url, "%{byte:02X}").unwrap();
        }
    }
    url
}

fn verify_existing(path: &Path, asset: &Asset) -> Result<()> {
    let metadata = fs::metadata(path)?;
    if !metadata.is_file() || metadata.len() != asset.size_bytes {
        bail!("Expected a file of {} bytes", asset.size_bytes);
    }
    copy_and_verify(File::open(path)?, io::sink(), asset)
}

fn copy_and_verify(reader: impl Read, mut writer: impl Write, asset: &Asset) -> Result<()> {
    let mut reader = reader.take(asset.size_bytes.saturating_add(1));
    let mut hasher = digest::Context::new(&SHA256);
    let mut buffer = [0; 64 * 1024];
    let mut bytes = 0;
    loop {
        let read = reader
            .read(&mut buffer)
            .context("Failed to read index bytes")?;
        if read == 0 {
            break;
        }
        let chunk = &buffer[..read];
        writer
            .write_all(chunk)
            .context("Failed to write index bytes")?;
        hasher.update(chunk);
        bytes += read as u64;
    }
    if bytes != asset.size_bytes {
        bail!(
            "Size mismatch: expected {} bytes, received {bytes}",
            asset.size_bytes
        );
    }
    let mut actual = String::with_capacity(64);
    for byte in hasher.finish().as_ref() {
        write!(actual, "{byte:02x}").unwrap();
    }
    if !actual.eq_ignore_ascii_case(&asset.sha256) {
        bail!(
            "SHA-256 mismatch: expected {}, computed {actual}",
            asset.sha256
        );
    }
    Ok(())
}

/// Treat existing directories and paths ending in a separator as directories
fn resolve_output_path(output: Option<&Path>, filename: &str) -> PathBuf {
    match output {
        Some(p) if p.is_dir() || p.to_string_lossy().ends_with(std::path::is_separator) => {
            p.join(filename)
        }
        Some(p) => p.to_path_buf(),
        None => PathBuf::from(filename),
    }
}
