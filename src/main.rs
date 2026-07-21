use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
#[cfg(feature = "fetch")]
use deacon::index_fetch;
use deacon::{
    index_diff, index_dump, index_filter, index_freeze, index_info, index_intersect, index_union,
    ComplexityAlgorithm, FilterConfig, IndexConfig, DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE,
};
use serde::{Deserialize, Serialize};
use std::io::{Read, Write};
use std::os::unix::net::{UnixListener, UnixStream};
use std::path::PathBuf;

#[derive(Parser, Serialize, Deserialize)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
    #[arg(long)]
    /// Execute command using an existing server process
    use_server: bool,
}

#[derive(Subcommand, Serialize, Deserialize)]
enum Command {
    /// Build, inspect, compose and fetch minimizer indexes
    Index {
        #[command(subcommand)]
        command: IndexCommand,
    },
    /// Retain or deplete sequence records with sufficient minimizer hits to an indexed query
    Filter {
        /// Path to minimizer index file
        index: PathBuf,

        /// Optional path to fastx file (or - for stdin)
        #[arg(default_value = "-")]
        input: String,

        /// Optional path to second paired fastx file
        input2: Option<String>,

        /// Minimum absolute number of minimizer hits for a match
        #[arg(short = 'a', long = "abs-threshold", default_value_t = 2, value_parser = clap::value_parser!(u16).range(1..))]
        abs_threshold: u16,

        /// Minimum relative proportion (0.0-1.0) of minimizer hits for a match
        #[arg(short = 'r', long = "rel-threshold", default_value_t = 0.01)]
        rel_threshold: f64,

        /// Search only the first N nucleotides per sequence (0 = entire sequence)
        #[arg(short = 'p', long = "prefix-length", default_value_t = 0)]
        prefix_length: usize,

        /// Ignore minimizers below this kdust complexity threshold (0.0-1.0)
        #[arg(short = 'c', long = "complexity-threshold")]
        complexity_threshold: Option<f32>,

        /// Discard matching sequences (invert filtering behaviour)
        #[arg(short = 'd', long = "deplete", default_value_t = false)]
        deplete: bool,

        /// Replace sequence headers with incrementing numbers
        #[arg(short = 'R', long = "rename", default_value_t = false)]
        rename: bool,

        /// Replace sequence headers with incrementing numbers and random suffixes
        #[arg(long = "rename-random", default_value_t = false)]
        rename_random: bool,

        /// Output FASTA format regardless of input format
        #[arg(short = 'f', long = "fasta", default_value_t = false)]
        output_fasta: bool,

        /// Path to output fastx file (stdout if not specified; detects .gz and .zst)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,

        /// Optional path to second paired output fastx file (detects .gz and .zst)
        #[arg(short = 'O', long = "output2")]
        output2: Option<String>,

        /// Path to JSON summary output file
        #[arg(short = 's', long = "summary")]
        summary: Option<PathBuf>,

        /// Number of threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: u16,

        /// Number of threads used for output compression (0 = auto)
        #[arg(long = "compression-threads", default_value_t = 0)]
        compression_threads: u16,

        /// Output compression level (1-9 for gz & xz; 1-22 for zstd)
        #[arg(long = "compression-level", default_value_t = 2)]
        compression_level: u8,

        /// Output sequences with minimizer hits to stderr
        #[arg(long = "debug", default_value_t = false)]
        debug: bool,

        /// Treat INPUT as interleaved paired reads from a file or stdin
        #[arg(
            long = "interleaved",
            default_value_t = false,
            conflicts_with = "input2"
        )]
        interleaved: bool,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,
    },
    /// Start/stop a server process for reduced latency filtering
    Server {
        #[command(subcommand)]
        command: ServerCommand,
    },
    /// Show citation information
    Cite,
}

#[derive(Subcommand, Serialize, Deserialize)]
enum ServerCommand {
    /// Start the server
    Start {
        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: u16,
    },
    /// Print whether a server is running and the current index.
    Status,
    /// Stop the running server
    Stop,
}

#[derive(Subcommand, Serialize, Deserialize)]
enum IndexCommand {
    /// Index minimizers contained within a fastx file
    Build {
        /// Path to input fastx file (or - for stdin; supports gz, zst and xz compression)
        input: PathBuf,

        /// K-mer length used for indexing (k+w-1 must be <= 96 and odd)
        #[arg(short = 'k', default_value_t = DEFAULT_KMER_LENGTH)]
        kmer_length: u8,

        /// Minimizer window size used for indexing
        #[arg(short = 'w', default_value_t = DEFAULT_WINDOW_SIZE)]
        window_size: u8,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: u16,

        /// Suppress sequence header output
        #[arg(short = 'q', long = "quiet")]
        quiet: bool,
    },
    /// Combine multiple minimizer indexes (A ∪ B…)
    Union {
        /// Path(s) to one or more index file(s)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
    /// Intersect multiple minimizer indexes (A ∩ B…)
    Intersect {
        /// Path(s) to two or more index file(s)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
    /// Subtract minimizers in one index from another (A - B)
    Diff {
        /// Path to first index file
        #[arg(required = true)]
        first: PathBuf,

        /// Path to second index file or FASTX file (or - for stdin when using FASTX)
        #[arg(required = true)]
        second: PathBuf,

        /// K-mer length (required if second argument is FASTX file, 1-32)
        #[arg(short = 'k', long = "kmer-length", value_parser = clap::value_parser!(u8).range(1..=32))]
        kmer_length: Option<u8>,

        /// Window size (required if second argument is FASTX file)
        #[arg(short = 'w', long = "window-size")]
        window_size: Option<u8>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: u16,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
    /// Dump minimizer index to fasta
    Dump {
        /// Path to index file
        index: PathBuf,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
    /// Discard minimizers below a complexity threshold
    Filter {
        /// Path to index file
        index: PathBuf,

        /// Complexity measure (~0.9 recommended for kdust)
        #[arg(short = 'a', long = "algorithm", value_enum, default_value_t = ComplexityAlgorithm::Kdust)]
        algorithm: ComplexityAlgorithm,

        /// Discard minimizers with complexity below this threshold (0.0-1.0)
        #[arg(short = 'c', long = "complexity-threshold")]
        threshold: f32,

        /// Invert: keep only minimizers below the threshold
        #[arg(short = 'i', long = "invert")]
        invert: bool,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
    /// Show index information
    Info {
        /// Path to index file
        index: PathBuf,
    },
    /// Freeze an index into a binary fuse filter (BFF) index (k<=32)
    Freeze {
        /// Path to input index file
        index: PathBuf,

        /// Path to output file (stdout if not specified)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,

        /// Fingerprint width in bits: 32 (FP ~2^-32, ~36 bits/key) or 16 (FP ~2^-16, ~18 bits/key)
        #[arg(short = 'b', long = "bits", default_value_t = 32, value_parser = parse_fingerprint_bits)]
        bits: u8,
    },
    /// Fetch a pre-built index from remote storage
    #[cfg(feature = "fetch")]
    Fetch {
        /// Index name (e.g., panhuman-1)
        #[arg(default_value = "panhuman-1")]
        index_name: String,

        /// K-mer length
        #[arg(short = 'k', default_value_t = DEFAULT_KMER_LENGTH)]
        kmer_length: u8,

        /// Minimizer window size
        #[arg(short = 'w', default_value_t = DEFAULT_WINDOW_SIZE)]
        window_size: u8,

        /// Path to output file (default: ./)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },
}

/// client -> server
#[derive(Serialize, Deserialize)]
enum Reply {
    /// Reply for `server status`.
    IndexPath(Option<PathBuf>),
    Done,
}

/// Parse and validate the BFF fingerprint width (16 or 32 bits)
fn parse_fingerprint_bits(s: &str) -> Result<u8, String> {
    match s {
        "16" => Ok(16),
        "32" => Ok(32),
        _ => Err(format!("expected 16 or 32, got `{s}`")),
    }
}

fn print_citation() {
    println!("Bede Constantinides, John Lees, Derrick W Crook.");
    println!("\"Deacon: fast sequence filtering and contaminant depletion\"");
    println!("bioRxiv 2025.06.09.658732");
    println!("https://doi.org/10.1101/2025.06.09.658732");
}

fn main() -> Result<()> {
    // Check we have either AVX2 or NEON
    #[cfg(not(any(target_feature = "avx2", target_feature = "neon")))]
    {
        eprintln!(
            "Warning: SIMD acceleration is unavailable. For best performance, compile with `cargo build --release -C target-cpu=native`"
        );
    }

    // If the binary was compiled with AVX2, check that the machine supports it at runtime.
    ensure_simd::ensure_simd();

    let cli = Cli::parse();

    // Start the server if requested.
    if let Command::Server {
        command: ServerCommand::Start { threads },
    } = &cli.command
    {
        // First try connecting to an existing server.
        let connect = UnixStream::connect("deacon_server_socket");
        if connect.is_ok() {
            return Err(anyhow::anyhow!("Server is already running."));
        }

        rayon::ThreadPoolBuilder::new()
            .num_threads(*threads as usize)
            .build_global()
            .context("Failed to initialize thread pool")?;

        // Remove existing socket if present
        let _ = std::fs::remove_file("deacon_server_socket");
        let listener = UnixListener::bind("deacon_server_socket")?;
        'stream: for stream in listener.incoming() {
            let mut stream = stream.unwrap();
            let mut message = vec![];
            let mut buf = vec![0; 10000];
            loop {
                let len = stream.read(&mut buf)?;
                if len == 0 {
                    // drop this message
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
            let message: Command = serde_json::from_slice(&message).unwrap();
            match message {
                Command::Server {
                    command: ServerCommand::Start { .. },
                } => {
                    // just reply Done from already-started server.
                    serde_json::to_writer(stream, &Reply::Done)?;
                }
                Command::Server {
                    command: ServerCommand::Status,
                } => {
                    let reply = Reply::IndexPath(deacon::current_index_path());
                    serde_json::to_writer(stream, &reply)?;
                }
                Command::Server {
                    command: ServerCommand::Stop,
                } => {
                    eprintln!("stopping the server");
                    serde_json::to_writer(stream, &Reply::Done)?;
                    let _ = std::fs::remove_file("deacon_server_socket");
                    break;
                }
                command => {
                    process_command(&command)?;
                    serde_json::to_writer(stream, &Reply::Done)?;
                }
            }
        }

        return Ok(());
    }

    // Send a command to the server.
    if cli.use_server || matches!(cli.command, Command::Server { .. }) {
        let mut stream = UnixStream::connect("deacon_server_socket")?;
        serde_json::to_writer(&stream, &cli.command)?;
        stream.write_all(b"\0")?;
        stream.flush()?;
        let message: Reply = serde_json::from_reader(stream).unwrap();
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
        }

        return Ok(());
    }

    process_command(&cli.command)?;

    Ok(())
}

fn process_command(command: &Command) -> Result<(), anyhow::Error> {
    match &command {
        Command::Cite => {
            print_citation();
        }
        Command::Server { .. } => {
            unreachable!("Server commands are handled before this function is called")
        }
        Command::Index { command } => match command {
            IndexCommand::Build {
                input,
                kmer_length,
                window_size,
                output,
                threads,
                quiet,
            } => {
                let config = IndexConfig {
                    input_path: input.clone(),
                    kmer_length: *kmer_length,
                    window_size: *window_size,
                    output_path: output.clone(),
                    threads: *threads,
                    quiet: *quiet,
                };
                config
                    .execute()
                    .context("Failed to run index build command")?;
            }
            IndexCommand::Info { index } => {
                index_info(index).context("Failed to run index info command")?;
            }
            #[cfg(feature = "fetch")]
            IndexCommand::Fetch {
                index_name,
                kmer_length,
                window_size,
                output,
            } => {
                index_fetch(index_name, *kmer_length, *window_size, output.as_deref())
                    .context("Failed to run index fetch command")?;
            }
            IndexCommand::Union { inputs, output } => {
                index_union(inputs, output.as_deref())
                    .context("Failed to run index union command")?;
            }
            IndexCommand::Intersect { inputs, output } => {
                index_intersect(inputs, output.as_deref())
                    .context("Failed to run index intersect command")?;
            }
            IndexCommand::Diff {
                first,
                second,
                kmer_length,
                window_size,
                output,
                threads,
            } => {
                index_diff(
                    first,
                    second,
                    *kmer_length,
                    *window_size,
                    *threads,
                    output.as_deref(),
                )
                .context("Failed to run index diff command")?;
            }
            IndexCommand::Dump { index, output } => {
                index_dump(index, output.as_deref()).context("Failed to run index dump command")?;
            }
            IndexCommand::Filter {
                index,
                algorithm,
                threshold,
                invert,
                output,
            } => {
                index_filter(index, output.as_deref(), *algorithm, *threshold, *invert)
                    .context("Failed to run index filter command")?;
            }
            IndexCommand::Freeze {
                index,
                output,
                bits,
            } => {
                index_freeze(index, output.as_deref(), *bits)
                    .context("Failed to run index freeze command")?;
            }
        },
        Command::Filter {
            index: minimizers,
            input,
            input2,
            interleaved,
            output,
            output2,
            abs_threshold,
            rel_threshold,
            prefix_length,
            complexity_threshold,
            summary,
            deplete,
            rename,
            rename_random,
            output_fasta,
            threads,
            compression_level,
            compression_threads,
            debug,
            quiet,
        } => {
            // Validate output2 usage
            if output2.is_some() && input2.is_none() && !interleaved {
                eprintln!(
                    "Warning: --output2 specified but no second input file provided. --output2 will be ignored."
                );
            }

            let config = FilterConfig {
                minimizers_path: minimizers,
                input_path: input,
                input2_path: input2.as_deref(),
                interleaved: *interleaved,
                output_path: output.as_ref().map(|p| p.as_path()),
                output2_path: output2.as_deref(),
                abs_threshold: *abs_threshold as usize,
                rel_threshold: *rel_threshold,
                prefix_length: *prefix_length,
                complexity_threshold: *complexity_threshold,
                summary_path: summary.as_ref(),
                deplete: *deplete,
                rename: *rename,
                rename_random: *rename_random,
                output_fasta: *output_fasta,
                threads: *threads,
                compression_level: *compression_level,
                compression_threads: *compression_threads,
                debug: *debug,
                quiet: *quiet,
            };
            config.execute().context("Failed to run filter command")?;
        }
    }

    Ok(())
}
