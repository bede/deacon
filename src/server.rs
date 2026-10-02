//! Binary server retaining an index between requests.

use crate::{Command, FilterArgs, ServerCommand, process_command};
use anyhow::{Context, Result};
use deacon::{FilterConfig, Index, filter_files, load_filter_index};
use serde::{Deserialize, Serialize};
use std::io::{Read, Write};
use std::os::unix::net::{UnixListener, UnixStream};
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;
use tracing::{info, warn};

const SOCKET: &str = "deacon_server_socket";

#[derive(Serialize, Deserialize)]
enum Reply {
    /// Reply for `server status`.
    IndexPath(Option<PathBuf>),
    Done,
    Error(String),
}

struct LoadedIndex {
    path: PathBuf,
    complexity_threshold: Option<f32>,
    index: Arc<Index>,
}

/// Filter with the cached index, reloading if its path changes.
fn filter(index: &mut Option<LoadedIndex>, args: &FilterArgs) -> Result<()> {
    // `server start` controls logging, so ignore the request's -q.
    let config = FilterConfig {
        progress: true,
        ..args.to_config()
    };
    let mut load_time = None;
    let loaded = index.as_ref().is_some_and(|i| {
        i.path == args.index && i.complexity_threshold == args.complexity_threshold
    });
    if !loaded {
        // Free the old index before loading.
        *index = None;
        let start = Instant::now();
        let loaded_index = load_filter_index(&args.index, args.complexity_threshold)?;
        load_time = Some(start.elapsed());
        *index = Some(LoadedIndex {
            path: args.index.clone(),
            complexity_threshold: args.complexity_threshold,
            index: Arc::new(loaded_index),
        });
    }
    let index = index.as_ref().unwrap();
    let label = args.index.to_string_lossy();
    filter_files(Arc::clone(&index.index), &label, load_time, &config)?;
    Ok(())
}

/// Serve until stopped.
pub fn start(threads: u16) -> Result<()> {
    // First try connecting to an existing server.
    if UnixStream::connect(SOCKET).is_ok() {
        return Err(anyhow::anyhow!("Server is already running."));
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads as usize)
        .build_global()
        .context("Failed to initialize thread pool")?;

    // Remove existing socket if present
    let _ = std::fs::remove_file(SOCKET);
    let listener = UnixListener::bind(SOCKET)?;
    let mut index: Option<LoadedIndex> = None;

    // Loop over incoming connections.
    'stream: for stream in listener.incoming() {
        let mut stream = match stream {
            Ok(s) => s,
            Err(e) => {
                warn!("Failed to accept incoming connection: {e}");
                continue 'stream;
            }
        };
        let mut message = vec![];
        let mut buf = vec![0; 10000];
        loop {
            let len = match stream.read(&mut buf) {
                Ok(len) => len,
                Err(e) => {
                    warn!("Failed to read request from client: {e}");
                    continue 'stream;
                }
            };
            if len == 0 {
                // drop this message
                warn!("Incoming request was empty");
                continue 'stream;
            }
            let buf = &buf[..len];
            message.extend_from_slice(buf);
            if buf.contains(&0) {
                assert_eq!(buf.last(), Some(&0));
                message.pop();
                break;
            }
        }
        let message: Command = match serde_json::from_slice(&message) {
            Ok(message) => message,
            Err(e) => {
                warn!("Failed to parse request from client: {e}");
                continue 'stream;
            }
        };
        let reply_status = match message {
            Command::Server {
                command: ServerCommand::Start { .. },
            } => {
                // just reply Done from already-started server.
                serde_json::to_writer(stream, &Reply::Done)
            }
            Command::Server {
                command: ServerCommand::Status,
            } => {
                let path = index.as_ref().map(|i| i.path.clone());
                serde_json::to_writer(stream, &Reply::IndexPath(path))
            }
            Command::Server {
                command: ServerCommand::Stop,
            } => {
                info!("Stopping the server");
                serde_json::to_writer(stream, &Reply::Done)?;
                let _ = std::fs::remove_file(SOCKET);
                break;
            }
            command => {
                let result = match command {
                    Command::Filter(args) => {
                        filter(&mut index, &args).context("Failed to run filter command")
                    }
                    command => process_command(command),
                };
                let reply = match result {
                    Ok(()) => Reply::Done,
                    Err(e) => Reply::Error(format!("{e:#}")),
                };
                serde_json::to_writer(stream, &reply)
            }
        };
        if let Err(e) = reply_status {
            warn!("Failed to send reply to client: {e}");
        }
    }

    Ok(())
}

/// Send a command and report the reply.
pub fn send(command: &Command) -> Result<()> {
    let mut stream = UnixStream::connect(SOCKET)?;
    serde_json::to_writer(&stream, command)?;
    stream.write_all(b"\0")?;
    stream.flush()?;
    let message: Reply = serde_json::from_reader(stream)
        .map_err(|e| anyhow::anyhow!("Could not read the server response:\n{e}"))?;
    match message {
        Reply::IndexPath(index_path) => {
            println!("Server is running.");
            if let Some(index_path) = index_path {
                println!("Current index: {}", index_path.display());
            } else {
                println!("No index is loaded yet.");
            }
        }
        Reply::Done => {}
        Reply::Error(e) => {
            return Err(anyhow::anyhow!(
                "The server had an error while processing the command:\n{e}"
            ));
        }
    }

    Ok(())
}
