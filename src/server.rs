//! Functionality to create a server endpoint which can be used to filter based on a pre-loaded index
use std::path::PathBuf;
use std::sync::{Mutex, OnceLock};

use crate::index::{IndexHeader, load_minimizer_hashes};
use crate::remote_filter::{paired_should_keep, unpaired_should_keep};
use crate::server_common::{FilterResponse, PairedFilterRequest, UnpairedFilterRequest};
use axum::{
    Json, Router,
    extract::DefaultBodyLimit,
    routing::{get, post},
};
use rustc_hash::FxHashSet;

/// Shared index file between endpoint calls.
/// Annoyingly, we have to use an Option as the default/empty FxHashSet is not static
static INDEX: Mutex<Option<FxHashSet<u64>>> = Mutex::new(None);

/// Shared index header between endpoint calls.
/// Initalised to a dummy value, which will be replaced when the index is loaded.
static INDEX_HEADER: Mutex<IndexHeader> = Mutex::new(IndexHeader {
    format_version: 0,
    kmer_length: 0,
    window_size: 0,
});

/// Shared index hash between endpoint calls.
static INDEX_HASH: Mutex<Option<String>> = Mutex::new(None);

/// Just for ensuring we get a single tracing setup.
/// Mostly needed as tests otherwise try to spawn multiple
static TRACING: OnceLock<()> = OnceLock::new();

/// Starts the server with the given index path and port.
/// To log the server's connections, set `RUST_LOG=trace` in your environment variables.
pub async fn run_server(index_path: PathBuf, port: u16) {
    // initialize tracing
    TRACING.get_or_init(|| {
        tracing_subscriber::fmt::init();
    });

    eprintln!("Loading index from: {}", index_path.display());
    // Load the index before starting the server to ensure it's available for requests
    load_index(index_path);
    eprintln!("Loaded index!");

    // build our application with a route
    let app = Router::new()
        // `GET /` goes to `root`
        .route("/", get(root))
        // `GET /index_header` returns the index header
        .route("/index_header", get(index_header))
        // `GET /index_version` returns the index version (hash)
        .route("/index_version", get(index_version))
        .route("/should_output_paired", post(should_output_paired))
        .route("/should_output_unpaired", post(should_output_unpaired))
        // Increase the body limit to 2GB to ensure we don't error on large payloads
        .layer(DefaultBodyLimit::max(2147483648));

    // run our app with hyper, listening globally
    let listener = tokio::net::TcpListener::bind("0.0.0.0:".to_owned() + &port.to_string())
        .await
        .unwrap();
    axum::serve(listener, app).await.unwrap();
}

/// Load the index from the specified path.
fn load_index(index_path: PathBuf) {
    // Load the hash as well as the file contents for returning as an ugly (but reliable) version
    let bytes = std::fs::read(index_path.clone()).unwrap();
    let hash = sha256::digest(&bytes);
    *INDEX_HASH.lock().unwrap() =
        Some(index_path.clone().into_os_string().into_string().unwrap() + "@" + &hash);

    let result = load_minimizer_hashes(&Some(&index_path), &None);
    match result {
        Ok((minimizers, header)) => {
            *INDEX.lock().unwrap() = minimizers;
            *INDEX_HEADER.lock().unwrap() = header;
        }
        Err(e) => {
            eprintln!("Failed to load index: {e}");
            std::process::exit(1);
        }
    }
}

/// Basic root, returing a message indicating the index is loaded
/// Endpoint is `/`
pub async fn root() -> String {
    let index = INDEX.lock();
    match index {
        Ok(index) => {
            let index = index.as_ref().expect("Index not loaded");
            let header = INDEX_HEADER.lock().unwrap();
            format!(
                "Index loaded with {} minimizers and header: {:?}",
                index.len(),
                header
            )
        }
        Err(e) => format!("Error accessing index: {e}"),
    }
}

/// Endpoint to return the header of the loaded index
/// Endpoint is `/index_header`
pub async fn index_header() -> Json<IndexHeader> {
    let header = INDEX_HEADER.lock().unwrap();
    Json(header.clone())
}

/// Endpoint to return the loaded index version
/// Endpoint is `/index_version`
pub async fn index_version() -> String {
    let hash = INDEX_HASH.lock().unwrap();
    hash.clone().unwrap()
}

async fn should_output_paired(Json(request): Json<PairedFilterRequest>) -> Json<FilterResponse> {
    // Quickly wrangle the seqs into slices from vecs as serde can't do it directly
    let input_minimizers_and_positions: Vec<(Vec<u64>, Vec<u32>, Vec<&[u8]>)> = request
        .input
        .iter()
        .map(|(minimizers, positions, seqs)| {
            (
                minimizers.to_vec(),
                positions.to_vec(),
                seqs.iter().map(|s| s.as_slice()).collect(),
            )
        })
        .collect();
    let index = INDEX.lock();
    match index {
        Ok(index) => {
            let index = index.as_ref().expect("Index not loaded");
            let result = paired_should_keep(
                &input_minimizers_and_positions,
                request.kmer_length,
                index,
                request.abs_threshold,
                request.rel_threshold,
                request.deplete,
                request.debug,
            );
            Json(FilterResponse {
                should_output: result,
            })
        }
        Err(e) => panic!("Error accessing index: {e}"),
    }
}

async fn should_output_unpaired(
    Json(request): Json<UnpairedFilterRequest>,
) -> Json<FilterResponse> {
    let index = INDEX.lock();
    match index {
        Ok(index) => {
            let index = index.as_ref().expect("Index not loaded");
            let result = unpaired_should_keep(
                &request.input,
                request.kmer_length,
                index,
                request.abs_threshold,
                request.rel_threshold,
                request.deplete,
                request.debug,
            );
            Json(FilterResponse {
                should_output: result,
            })
        }
        Err(e) => panic!("Error accessing index: {e}"),
    }
}
