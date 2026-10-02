use crate::{Command, ServerCommand, process_command};
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::io::{Read, Write};
use std::os::unix::net::{UnixListener, UnixStream};
use std::path::PathBuf;
use tracing::{info, warn};

/// client -> server
#[derive(Serialize, Deserialize)]
enum Reply {
    /// Reply for `server status`.
    IndexPath(Option<PathBuf>),
    Done,
    Error(String),
}

pub fn start(threads: u16) -> Result<()> {
    // First try connecting to an existing server.
    let connect = UnixStream::connect("deacon_server_socket");
    if connect.is_ok() {
        return Err(anyhow::anyhow!("Server is already running."));
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads as usize)
        .build_global()
        .context("Failed to initialize thread pool")?;

    // Remove existing socket if present
    let _ = std::fs::remove_file("deacon_server_socket");
    let listener = UnixListener::bind("deacon_server_socket")?;

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
                let reply = Reply::IndexPath(deacon::current_index_path());
                serde_json::to_writer(stream, &reply)
            }
            Command::Server {
                command: ServerCommand::Stop,
            } => {
                info!("Stopping the server");
                serde_json::to_writer(stream, &Reply::Done)?;
                let _ = std::fs::remove_file("deacon_server_socket");
                break;
            }
            command => {
                let result = process_command(command);
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

pub fn send(command: &Command) -> Result<()> {
    let mut stream = UnixStream::connect("deacon_server_socket")?;
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
