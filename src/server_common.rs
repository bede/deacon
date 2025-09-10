//! Common structures and types used in the client and server
use crate::IndexHeader;
use anyhow::Result;
use reqwest::blocking::Client;
use serde::{Deserialize, Serialize};

/// Request structure for filtering minimizers unpaired reads
#[derive(Serialize, Deserialize)]
pub struct UnpairedFilterRequest {
    /// Prehashed minimizers for input
    /// Tuple of (minimizer hashes, positions, effective sequences)
    pub input: Vec<(Vec<u64>, Vec<u32>, Vec<u8>)>,

    /// Mininum number (integer) of minimizer hits for a match
    pub abs_threshold: usize,

    /// Mininum proportion (float) of minimizer hits for a match
    pub rel_threshold: f64,

    /// Whether running in deplete mode
    pub deplete: bool,

    /// kmer length used to compute minimizers
    pub kmer_length: u8,

    /// Whether running in debug mode
    pub debug: bool,
}

/// Request structure for filtering minimizers from paired reads
#[derive(Serialize, Deserialize)]
pub struct PairedFilterRequest {
    /// Prehashed minimizers for input
    /// Tuple of (minimizer hashes, positions, effective sequences)
    pub input: Vec<(Vec<u64>, Vec<u32>, Vec<Vec<u8>>)>,

    /// Mininum number (integer) of minimizer hits for a match
    pub abs_threshold: usize,

    /// Mininum proportion (float) of minimizer hits for a match
    pub rel_threshold: f64,

    /// Whether running in deplete mode
    pub deplete: bool,

    /// kmer length used to compute minimizers
    pub kmer_length: u8,

    /// Whether running in debug mode
    pub debug: bool,
}

/// Response structure for filter results
/// Returns whether this set of minimizers should be output
#[derive(Serialize, Deserialize)]
pub struct FilterResponse {
    /// Indicates whether this set of minimizers should be output
    /// Tuple of (should_keep, hit_count, total_minimizers, hit_kmers)
    pub should_output: Vec<(bool, usize, usize, Vec<String>)>,
}

/// Get the header of the index loaded into a remote server
/// Required in order to ensure that the locally computed minimizers match
/// the kmer length and window size
pub fn get_server_index_header(server_address: &str) -> Result<IndexHeader> {
    // Create a client to send the minimizers to the server
    let client = Client::new();

    // Send the minimizers as a POST request
    let response = client
        .get(server_address.to_owned() + "/index_header")
        .send()?;

    // Check if the response indicates a match
    if response.status().is_success() {
        Ok(response.json::<IndexHeader>()?)
    } else {
        Err(anyhow::anyhow!(
            "Server returned an error: {}",
            response.status()
        ))
    }
}
