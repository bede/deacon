use crate::minimizers::validate_k_w;
use anyhow::{Context, Result, bail};
use serde::Deserialize;
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;

#[derive(Deserialize)]
pub(super) struct Manifest {
    pub(super) default_index: Option<String>,
    indexes: BTreeMap<String, IndexGroup>,
    #[serde(default)]
    aliases: BTreeMap<String, String>,
}

#[derive(Deserialize)]
struct IndexGroup {
    k: u8,
    w: u8,
    default_asset: Option<String>,
    assets: BTreeMap<String, Asset>,
}

#[derive(Deserialize)]
pub(super) struct Asset {
    pub(super) size_bytes: u64,
    pub(super) sha256: String,
    #[serde(default)]
    opts: BTreeMap<String, Value>,
}

impl IndexGroup {
    fn default_filename(&self) -> Option<&str> {
        if let Some(filename) = &self.default_asset {
            return Some(filename);
        }
        if self.assets.len() == 1 {
            return self.assets.keys().next().map(String::as_str);
        }
        None
    }
}

impl Manifest {
    pub(super) fn parse(body: &[u8]) -> Result<Self> {
        let document: Value =
            serde_json::from_slice(body).context("Invalid JSON index manifest")?;
        let version = document
            .get("schema_version")
            .and_then(Value::as_u64)
            .context("Expected an integer schema_version in the index manifest")?;
        if version != 1 {
            bail!("Unsupported index manifest schema version {version}");
        }
        let manifest: Self = serde_json::from_value(document).context("Invalid index manifest")?;
        manifest.validate()?;
        Ok(manifest)
    }

    fn validate(&self) -> Result<()> {
        let mut names = BTreeSet::new();
        for (name, group) in &self.indexes {
            validate_name(name)?;
            names.insert(name);
            validate_k_w(group.k, group.w)
                .with_context(|| format!("Invalid index group {name:?}"))?;
            let default = group.default_filename().with_context(|| {
                format!("Index group {name:?} needs one asset or an explicit default_asset")
            })?;
            if !group.assets.contains_key(default) {
                bail!("Default asset {default:?} is missing from index group {name:?}");
            }
        }
        for group in self.indexes.values() {
            for (filename, asset) in &group.assets {
                validate_filename(filename)?;
                if !names.insert(filename) {
                    bail!("Duplicate index name {filename:?}");
                }
                if asset.size_bytes == 0 {
                    bail!("Invalid size for {filename:?}: expected a nonzero size");
                }
                if asset.sha256.len() != 64
                    || !asset.sha256.bytes().all(|byte| byte.is_ascii_hexdigit())
                {
                    bail!("Invalid SHA-256 for {filename:?}");
                }
                if asset.opts.contains_key("k") || asset.opts.contains_key("w") {
                    bail!(
                        "Invalid options for {filename:?}: k and w are reserved index group fields"
                    );
                }
            }
        }
        for (alias, name) in &self.aliases {
            validate_name(alias)?;
            if !names.insert(alias) {
                bail!("Duplicate index name {alias:?}");
            }
            if !self.indexes.contains_key(name) {
                bail!("Alias {alias:?} points to unknown index group {name:?}");
            }
        }
        if let Some(name) = &self.default_index
            && self.find_asset(name).is_none()
        {
            bail!("Manifest default_index must name an asset filename, got {name:?}");
        }
        Ok(())
    }

    pub(super) fn resolve(&self, name: &str) -> Result<(&str, &Asset)> {
        let group_name = self.aliases.get(name).map_or(name, String::as_str);
        if let Some(group) = self.indexes.get(group_name) {
            let filename = group.default_filename().expect("validated default asset");
            return Ok((filename, &group.assets[filename]));
        }
        self.find_asset(name).with_context(|| {
            format!("Unknown index {name:?}. Run 'deacon index list' to see available names")
        })
    }

    fn find_asset(&self, filename: &str) -> Option<(&str, &Asset)> {
        self.indexes
            .values()
            .find_map(|group| group.assets.get_key_value(filename))
            .map(|(filename, asset)| (filename.as_str(), asset))
    }

    pub(super) fn write_listing(&self, mut output: impl Write) -> Result<()> {
        let mut groups: Vec<_> = self.indexes.iter().collect();
        let default_filename = self.default_index.as_deref();
        if let Some(filename) = default_filename {
            let position = groups
                .iter()
                .position(|(_, group)| group.assets.contains_key(filename))
                .expect("validated default asset belongs to a group");
            groups[..=position].rotate_right(1);
        }
        for (position, (name, group)) in groups.into_iter().enumerate() {
            let mut details = Vec::new();
            let aliases: Vec<_> = self
                .aliases
                .iter()
                .filter_map(|(alias, target)| (target == name).then_some(alias.as_str()))
                .collect();
            if !aliases.is_empty() {
                details.push(format!("aliases: {}", aliases.join(", ")));
            }
            if position == 0
                && let Some(filename) = default_filename
            {
                let variant = filename
                    .strip_prefix(name.as_str())
                    .filter(|suffix| matches!(*suffix, ".idx" | ".pidx"))
                    .unwrap_or(filename);
                details.push(format!("default fetch: {variant}"));
            }
            if details.is_empty() {
                writeln!(output, "{name}")?;
            } else {
                writeln!(output, "{name} ({})", details.join(", "))?;
            }
            for (filename, asset) in &group.assets {
                let size = format_size(asset.size_bytes);
                write!(
                    output,
                    "    {filename:<28}  {size:>10}  k={} w={}",
                    group.k, group.w
                )?;
                if let Some(value) = asset.opts.get("c") {
                    write!(output, " c={value}")?;
                }
                for (key, value) in asset.opts.iter().filter(|(key, _)| key.as_str() != "c") {
                    if !key.is_empty()
                        && key
                            .chars()
                            .all(|c| c.is_alphanumeric() || matches!(c, '_' | '-' | '.'))
                    {
                        write!(output, " {key}={value}")?;
                    } else {
                        write!(output, " {}={value}", serde_json::to_string(key)?)?;
                    }
                }
                writeln!(output)?;
            }
        }
        Ok(())
    }
}

fn format_size(bytes: u64) -> String {
    let units = ["B", "kB", "MB", "GB", "TB", "PB", "EB"];
    let mut size = bytes as f64;
    let mut unit = 0;
    while size >= 1000.0 {
        size /= 1000.0;
        unit += 1;
    }
    let precision = usize::from(unit != 0);
    format!("{size:.precision$} {}", units[unit])
}

fn validate_name(name: &str) -> Result<()> {
    if name.trim() != name || name.is_empty() || name.chars().any(char::is_control) {
        bail!("Invalid index name {name:?}");
    }
    Ok(())
}

pub(super) fn validate_filename(filename: &str) -> Result<()> {
    validate_name(filename)?;
    // Require a single portable path component
    if filename.ends_with('.') || filename.contains(['/', '\\', ':', '*', '?', '"', '<', '>', '|'])
    {
        bail!("Invalid index asset filename {filename:?}");
    }
    Ok(())
}
