//! Common structures and types used in the client and server
use crate::{IndexHeader, MatchThreshold};
use anyhow::Result;
use reqwest::blocking::Client;
use serde::{Deserialize, Serialize};

/// Request structure for filtering minimizers
#[derive(Serialize, Deserialize)]
pub struct FilterRequest {
    /// Prehashed minimizers for input
    pub input: Vec<Vec<u64>>,

    /// Mininum number (integer) or proportion (float) of minimizer hits for a match
    pub match_threshold: MatchThreshold,

    /// Whether running in deplete mode
    pub deplete: bool,
}

/// Response structure for filter results
/// Returns whether this set of minimizers should be output
#[derive(Serialize, Deserialize)]
pub struct FilterResponse {
    /// Indicates if the input minimizers should be output
    pub should_output: Vec<bool>,
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
