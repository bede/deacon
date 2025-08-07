use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use deacon::{
    DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE, MatchThreshold, build_index, diff_index, index_info,
    run_filter, union_index,
};
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Create and compose minimizer indexes
    Index {
        #[command(subcommand)]
        command: IndexCommands,
    },
    /// Keep or discard fastx records with sufficient minimizer hits to the index
    Filter {
        /// Path to minimizer index file
        index: PathBuf,

        /// Optional path to fastx file (or - for stdin)
        #[arg(default_value = "-")]
        input: String,

        /// Optional path to second paired fastx file (or - for interleaved stdin)
        input2: Option<String>,

        /// Path to output fastx file (or - for stdout; detects .gz and .zst)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Optional path to second paired output fastx file (detects .gz and .zst)
        #[arg(short = 'O', long = "output2")]
        output2: Option<String>,

        /// Mininum number (integer) or proportion (float) of minimizer hits for a match
        #[arg(short = 'm', long = "matches", default_value_t = MatchThreshold::Absolute(2))]
        match_threshold: MatchThreshold,

        /// Search only the first N nucleotides per sequence (0 = entire sequence)
        #[arg(short = 'p', long = "prefix-length", default_value_t = 0)]
        prefix_length: usize,

        /// Discard matching sequences (invert filtering behaviour)
        #[arg(short = 'd', long = "deplete", default_value_t = false)]
        deplete: bool,

        /// Replace sequence headers with incrementing numbers
        #[arg(short = 'r', long = "rename", default_value_t = false)]
        rename: bool,

        /// Path to JSON summary output file
        #[arg(short = 's', long = "summary")]
        summary: Option<PathBuf>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Output compression level (1-9 for gz & xz; 1-22 for zstd)
        #[arg(long = "compression-level", default_value_t = 2)]
        compression_level: u8,
    },

    Server {
        /// Path to minimizer index file
        index: PathBuf,

        /// Port to run the server on
        #[arg(short = 'p', long = "port", default_value_t = 8888)]
        port: u16,
    },

    /// Alternate version of Filter, swapping local compute for passing to a server
    /// which has the index pre-loaded. Will inevitably be slower than local filtering,
    /// but saves on index loading. Better used for cases of small input + large index
    Client {
        /// Server address to connect to (including port)
        server_address: String,

        /// Optional path to fastx file (or - for stdin)
        #[arg(default_value = "-")]
        input: String,

        /// Optional path to second paired fastx file (or - for interleaved stdin)
        input2: Option<String>,

        /// Path to output fastx file (or - for stdout; detects .gz and .zst)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Optional path to second paired output fastx file (detects .gz and .zst)
        #[arg(short = 'O', long = "output2")]
        output2: Option<String>,

        /// Mininum number (integer) or proportion (float) of minimizer hits for a match
        #[arg(short = 'm', long = "matches", default_value_t = MatchThreshold::Absolute(2))]
        match_threshold: MatchThreshold,

        /// Search only the first N nucleotides per sequence (0 = entire sequence)
        #[arg(short = 'p', long = "prefix-length", default_value_t = 0)]
        prefix_length: usize,

        /// Discard matching sequences (invert filtering behaviour)
        #[arg(short = 'd', long = "deplete", default_value_t = false)]
        deplete: bool,

        /// Replace sequence headers with incrementing numbers
        #[arg(short = 'r', long = "rename", default_value_t = false)]
        rename: bool,

        /// Path to JSON summary output file
        #[arg(short = 's', long = "summary")]
        summary: Option<PathBuf>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Output compression level (1-9 for gz & xz; 1-22 for zstd)
        #[arg(long = "compression-level", default_value_t = 2)]
        compression_level: u8,
    },
}

#[derive(Subcommand)]
enum IndexCommands {
    /// Index minimizers contained within a fastx file
    Build {
        /// Path to input fastx file (supports gz, zst and xz compression)
        input: PathBuf,

        /// K-mer length used for indexing
        #[arg(short = 'k', default_value_t = DEFAULT_KMER_LENGTH)]
        kmer_length: usize,

        /// Minimizer window size used for indexing
        #[arg(short = 'w', default_value_t = DEFAULT_WINDOW_SIZE)]
        window_size: usize,

        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Preallocated index capacity in millions of minimizers
        #[arg(short = 'c', long = "capacity", default_value_t = 400)]
        capacity_millions: usize,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,
    },
    /// Show index information
    Info {
        /// Path to index file
        index: PathBuf,
    },
    /// Combine multiple minimizer indexes (A ∪ B…)
    Union {
        /// Path(s) to one or more index file(s)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,

        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
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

        /// K-mer length (required if second argument is FASTX file)
        #[arg(short = 'k', long = "kmer-length")]
        kmer_length: Option<usize>,

        /// Window size (required if second argument is FASTX file)
        #[arg(short = 'w', long = "window-size")]
        window_size: Option<usize>,

        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: Option<PathBuf>,
    },
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

    match &cli.command {
        Commands::Index { command } => match command {
            IndexCommands::Build {
                input,
                kmer_length,
                window_size,
                output,
                capacity_millions,
                threads,
            } => {
                // Convert output string to Option<PathBuf>
                let output_path = if output == "-" {
                    None
                } else {
                    Some(PathBuf::from(output))
                };

                build_index(
                    input,
                    *kmer_length,
                    *window_size,
                    output_path,
                    *capacity_millions,
                    *threads,
                )
                .context("Failed to run index build command")?;
            }
            IndexCommands::Info { index } => {
                index_info(index).context("Failed to run index info command")?;
            }
            IndexCommands::Union { inputs, output } => {
                union_index(inputs, output.as_ref())
                    .context("Failed to run index union command")?;
            }
            IndexCommands::Diff {
                first,
                second,
                kmer_length,
                window_size,
                output,
            } => {
                diff_index(first, second, *kmer_length, *window_size, output.as_ref())
                    .context("Failed to run index diff command")?;
            }
        },
        Commands::Filter {
            index: minimizers,
            input,
            input2,
            output,
            output2,
            match_threshold,
            prefix_length,
            summary,
            deplete,
            rename,
            threads,
            compression_level,
        } => {
            // Validate output2 usage
            if output2.is_some() && input2.is_none() {
                eprintln!(
                    "Warning: --output2 specified but no second input file provided. --output2 will be ignored."
                );
            }

            run_filter(
                Some(minimizers),
                input,
                input2.as_deref(),
                output,
                output2.as_deref(),
                match_threshold,
                *prefix_length,
                summary.as_ref(),
                *deplete,
                *rename,
                *threads,
                *compression_level,
                None,
            )
            .context("Failed to run filter command")?;
        }
        Commands::Server { index, port } => {
            #[cfg(feature = "server")]
            {
                // Server needs to run async, so spawn an async runtime to run it
                let rt = tokio::runtime::Runtime::new().unwrap();
                rt.block_on(async {
                    eprintln!("Loading server!");
                    deacon::server::run_server(index.clone(), *port).await;
                });
            }
            #[cfg(not(feature = "server"))]
            {
                eprintln!(
                    "Server functionality is not enabled in this build. Please compile with the 'server' feature: `cargo build --features server`"
                );
                // Suppress dead code warning so this compiles without issue when server is not enabled
                let _ = (index, port);
                std::process::exit(1);
            }
        }
        Commands::Client {
            server_address,
            input,
            input2,
            output,
            output2,
            match_threshold,
            prefix_length,
            summary,
            deplete,
            rename,
            threads,
            compression_level,
        } => {
            #[cfg(feature = "server")]
            {
                // Validate output2 usage
                if output2.is_some() && input2.is_none() {
                    eprintln!(
                        "Warning: --output2 specified but no second input file provided. --output2 will be ignored."
                    );
                }
                let dummy_index: Option<&PathBuf> = None;
                run_filter(
                    dummy_index,
                    input,
                    input2.as_deref(),
                    output,
                    output2.as_deref(),
                    match_threshold,
                    *prefix_length,
                    summary.as_ref(),
                    *deplete,
                    *rename,
                    *threads,
                    *compression_level,
                    Some(server_address.to_string()),
                )
                .context("Failed to run filter with client functionality")?;
            }
            #[cfg(not(feature = "server"))]
            {
                eprintln!(
                    "Client functionality is not enabled in this build. Please compile with the 'server' feature: `cargo build --features server`"
                );
                // Suppress dead code warning so this compiles without issue when server is not enabled
                let _ = (
                    server_address,
                    input,
                    input2,
                    output,
                    output2,
                    match_threshold,
                    prefix_length,
                    summary,
                    deplete,
                    rename,
                    threads,
                    compression_level,
                );
                std::process::exit(1);
            }
        }
    }

    Ok(())
}
