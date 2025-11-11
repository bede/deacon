use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use deacon::{
    DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE, FilterConfig, IndexConfig, diff_index, dump_index,
    index_info, intersect_index, union_index,
};
use serde::{Deserialize, Serialize};
use std::io::{Read, Write};
use std::os::unix::net::{UnixListener, UnixStream};
use std::path::PathBuf;

#[cfg(feature = "fetch")]
use indicatif::ProgressBar;

#[derive(Parser, Serialize, Deserialize)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
    #[arg(long)]
    /// Execute command using an existing server process
    use_server: bool,
}

#[derive(Subcommand, Serialize, Deserialize)]
enum Commands {
    /// Build, compose and inspect minimizer indexes
    Index {
        #[command(subcommand)]
        command: IndexCommands,
    },
    /// Retain or deplete sequence records with sufficient minimizer hits to an indexed query
    Filter {
        /// Path to minimizer index file
        index: PathBuf,

        /// Optional path to fastx file (or - for stdin)
        #[arg(default_value = "-")]
        input: String,

        /// Optional path to second paired fastx file (or - for interleaved stdin)
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

        /// Discard matching sequences (invert filtering behaviour)
        #[arg(short = 'd', long = "deplete", default_value_t = false)]
        deplete: bool,

        /// Replace sequence headers with incrementing numbers
        #[arg(short = 'R', long = "rename", default_value_t = false)]
        rename: bool,

        /// Path to output fastx file (stdout if not specified; detects .gz and .zst)
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,

        /// Optional path to second paired output fastx file (detects .gz and .zst)
        #[arg(short = 'O', long = "output2")]
        output2: Option<String>,

        /// Path to JSON summary output file
        #[arg(short = 's', long = "summary")]
        summary: Option<PathBuf>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: u16,

        /// Number of threads for compression (0 = auto)
        #[arg(long = "compression-threads", default_value_t = 0)]
        compression_threads: u16,

        /// Output compression level (1-9 for gz & xz; 1-22 for zstd)
        #[arg(long = "compression-level", default_value_t = 2)]
        compression_level: u8,

        /// Output sequences with minimizer hits to stderr
        #[arg(long = "debug", default_value_t = false)]
        debug: bool,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,
    },
    /// Start/stop a server process for reduced latency filtering
    Server {
        #[command(subcommand)]
        command: ServerCommands,
    },
    /// Show citation information
    Cite,
}

#[derive(Subcommand, Serialize, Deserialize)]
enum ServerCommands {
    /// Start the server
    Start {
        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: u16,
    },
    /// Stop the running server
    Stop,
}

#[derive(Subcommand, Serialize, Deserialize)]
enum IndexCommands {
    /// Index minimizers contained within a fastx file
    Build {
        /// Path to input fastx file (supports gz, zst and xz compression)
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

        /// Minimum scaled entropy threshold for k-mer filtering (0.0-1.0)
        #[arg(short = 'e', long = "entropy-threshold", default_value = "0.0")]
        entropy_threshold: f32,
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
    /// Show index information
    Info {
        /// Path to index file
        index: PathBuf,
    },
    /// Fetch a pre-built index from remote storage
    #[cfg(feature = "fetch")]
    Fetch {
        /// Index identifier (e.g., panhuman-1)
        #[arg(default_value = "panhuman-1")]
        index_id: String,

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

#[derive(Serialize, Deserialize)]
enum Message {
    /// client -> server
    Command(Commands),
    /// server -> client
    Done,
}

fn print_citation() {
    println!("Bede Constantinides, John Lees, Derrick W Crook.");
    println!("\"Deacon: fast sequence filtering and contaminant depletion\"");
    println!("bioRxiv 2025.06.09.658732");
    println!("https://doi.org/10.1101/2025.06.09.658732");
}

#[cfg(feature = "fetch")]
fn fetch_index(
    index_id: &str,
    kmer_length: u8,
    window_size: u8,
    output: Option<&std::path::Path>,
) -> Result<()> {
    const BASE_URL: &str =
        "https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o";

    let filename = format!("{}.k{}w{}.idx", index_id, kmer_length, window_size);
    let url = format!("{}/deacon/{}", BASE_URL, filename);

    eprintln!("Fetching {}", url);

    let response = ureq::get(&url).call().context("Failed to download index")?;

    if response.status() != 200 {
        anyhow::bail!("Failed to fetch index: HTTP {}", response.status());
    }

    let content_length = response
        .headers()
        .get("Content-Length")
        .and_then(|v| v.to_str().ok())
        .and_then(|s| s.parse::<u64>().ok());

    let pb = ProgressBar::new(content_length.unwrap_or(0));

    let mut body = response.into_body();
    let output_path = output
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from(&filename));

    let mut temp_path = output_path.clone();
    temp_path.as_mut_os_string().push(".tmp");

    let mut file = std::fs::File::create(&temp_path).context("Failed to create temporary file")?;
    std::io::copy(&mut pb.wrap_read(body.as_reader()), &mut file)
        .context("Failed to write index to file")?;

    std::fs::rename(&temp_path, &output_path)
        .context("Failed to move downloaded file to final location")?;

    pb.finish_and_clear();
    eprintln!("Index saved to: {}", output_path.display());

    Ok(())
}

fn main() -> Result<()> {
    // Check we have either AVX2 or NEON
    #[cfg(not(any(target_feature = "avx2", target_feature = "neon")))]
    {
        eprintln!(
            "Warning: SIMD acceleration is unavailable. For best performance, compile with `cargo build --release -C target-cpu=native`"
        );
    }

    let cli = Cli::parse();

    if let Commands::Server { command } = &cli.command {
        if !cli.use_server {
            // Running server commands directly (not via --use-server)
            match command {
                ServerCommands::Start { threads } => {
                    rayon::ThreadPoolBuilder::new()
                        .num_threads(*threads as usize)
                        .build_global()
                        .context("Failed to initialize thread pool")?;

                    // Remove existing socket if present
                    let _ = std::fs::remove_file("deacon_server_socket");
                    let listener = UnixListener::bind("deacon_server_socket")?;
                    for stream in listener.incoming() {
                        let mut stream = stream.unwrap();
                        let mut message = vec![];
                        let mut buf = vec![0; 10000];
                        loop {
                            let len = stream.read(&mut buf)?;
                            let buf = &buf[..len];
                            message.extend_from_slice(buf);
                            if buf.contains(&0) {
                                assert_eq!(buf.last(), Some(&0));
                                message.pop();
                                break;
                            }
                        }
                        let message: Message = serde_json::from_slice(&message).unwrap();
                        match message {
                            Message::Command(Commands::Server {
                                command: ServerCommands::Stop,
                            }) => {
                                serde_json::to_writer(stream, &Message::Done)?;
                                let _ = std::fs::remove_file("deacon_server_socket");
                                break;
                            }
                            Message::Command(commands) => {
                                process_command(&commands)?;
                                serde_json::to_writer(stream, &Message::Done)?;
                            }
                            Message::Done => {
                                unreachable!("Server should not receive `Done` messages.")
                            }
                        }
                    }

                    return Ok(());
                }
                ServerCommands::Stop => {
                    panic!("Use `deacon --use-server server stop` to stop the server.")
                }
            }
        }
        // If --use-server is set, fall through to client code below
    }

    if cli.use_server {
        let mut stream = UnixStream::connect("deacon_server_socket")?;
        serde_json::to_writer(&stream, &Message::Command(cli.command))?;
        stream.write(b"\0")?;
        stream.flush()?;
        let message: Message = serde_json::from_reader(stream).unwrap();
        match message {
            Message::Done => {}
            _ => unreachable!("The client only expects to receive `Done` messages."),
        }

        return Ok(());
    }

    process_command(&cli.command)?;

    Ok(())
}

fn process_command(command: &Commands) -> Result<(), anyhow::Error> {
    match &command {
        Commands::Cite => {
            print_citation();
        }
        Commands::Server { .. } => {
            unreachable!("Server commands are handled before this function is called")
        }
        Commands::Index { command } => match command {
            IndexCommands::Build {
                input,
                kmer_length,
                window_size,
                output,
                threads,
                quiet,
                entropy_threshold,
            } => {
                let config = IndexConfig {
                    input_path: input.clone(),
                    kmer_length: *kmer_length,
                    window_size: *window_size,
                    output_path: output.clone(),
                    threads: *threads,
                    quiet: *quiet,
                    entropy_threshold: *entropy_threshold,
                };
                config
                    .execute()
                    .context("Failed to run index build command")?;
            }
            IndexCommands::Info { index } => {
                index_info(index).context("Failed to run index info command")?;
            }
            #[cfg(feature = "fetch")]
            IndexCommands::Fetch {
                index_id,
                kmer_length,
                window_size,
                output,
            } => {
                fetch_index(index_id, *kmer_length, *window_size, output.as_deref())
                    .context("Failed to run index fetch command")?;
            }
            IndexCommands::Union { inputs, output } => {
                union_index(inputs, output.as_deref())
                    .context("Failed to run index union command")?;
            }
            IndexCommands::Intersect { inputs, output } => {
                intersect_index(inputs, output.as_deref())
                    .context("Failed to run index intersect command")?;
            }
            IndexCommands::Diff {
                first,
                second,
                kmer_length,
                window_size,
                output,
                threads,
            } => {
                diff_index(
                    first,
                    second,
                    *kmer_length,
                    *window_size,
                    *threads,
                    output.as_deref(),
                )
                .context("Failed to run index diff command")?;
            }
            IndexCommands::Dump { index, output } => {
                dump_index(index, output.as_deref()).context("Failed to run index dump command")?;
            }
        },
        Commands::Filter {
            index: minimizers,
            input,
            input2,
            output,
            output2,
            abs_threshold,
            rel_threshold,
            prefix_length,
            summary,
            deplete,
            rename,
            threads,
            compression_level,
            compression_threads,
            debug,
            quiet,
        } => {
            // Validate output2 usage
            if output2.is_some() && input2.is_none() {
                eprintln!(
                    "Warning: --output2 specified but no second input file provided. --output2 will be ignored."
                );
            }

            let config = FilterConfig {
                minimizers_path: minimizers,
                input_path: input,
                input2_path: input2.as_deref(),
                output_path: output.as_ref().map(|p| p.as_path()),
                output2_path: output2.as_deref(),
                abs_threshold: *abs_threshold as usize,
                rel_threshold: *rel_threshold,
                prefix_length: *prefix_length,
                summary_path: summary.as_ref(),
                deplete: *deplete,
                rename: *rename,
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
