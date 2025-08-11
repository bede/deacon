use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use deacon::{
    DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE, build_index, diff_index, index_info, run_filter,
    union_index,
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
    /// Keep or discard DNA fastx records with sufficient minimizer hits to the index
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

        /// Path to JSON summary output file
        #[arg(short = 's', long = "summary")]
        summary: Option<PathBuf>,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

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
}

#[derive(Subcommand)]
enum IndexCommands {
    /// Index minimizers contained within a fastx file
    Build {
        /// Path to input fastx file (supports gz, zst and xz compression)
        input: PathBuf,

        /// K-mer length used for indexing (1-32)
        #[arg(short = 'k', default_value_t = DEFAULT_KMER_LENGTH, value_parser = clap::value_parser!(u8).range(1..=32))]
        kmer_length: u8,

        /// Minimizer window size used for indexing
        #[arg(short = 'w', default_value_t = DEFAULT_WINDOW_SIZE)]
        window_size: u8,

        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Preallocated index capacity in millions of minimizers
        #[arg(short = 'c', long = "capacity", default_value_t = 400)]
        capacity_millions: usize,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Suppress sequence header output
        #[arg(short = 'q', long = "quiet")]
        quiet: bool,

        /// Minimum scaled entropy threshold for k-mer filtering (0.0-1.0)
        #[arg(short = 'e', long = "entropy-threshold")]
        entropy_threshold: Option<f32>,
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

        /// Preallocated index capacity in millions of minimizers (overrides sum-based allocation)
        #[arg(short = 'c', long = "capacity")]
        capacity_millions: Option<usize>,
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
                quiet,
                entropy_threshold,
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
                    *quiet,
                    *entropy_threshold,
                )
                .context("Failed to run index build command")?;
            }
            IndexCommands::Info { index } => {
                index_info(index).context("Failed to run index info command")?;
            }
            IndexCommands::Union {
                inputs,
                output,
                capacity_millions,
            } => {
                union_index(inputs, output.as_ref(), *capacity_millions)
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
            abs_threshold,
            rel_threshold,
            prefix_length,
            summary,
            deplete,
            rename,
            threads,
            compression_level,
            debug,
            quiet,
        } => {
            // Validate output2 usage
            if output2.is_some() && input2.is_none() {
                eprintln!(
                    "Warning: --output2 specified but no second input file provided. --output2 will be ignored."
                );
            }

            run_filter(
                minimizers,
                &input,
                input2.as_deref(),
                output,
                output2.as_deref(),
                *abs_threshold as usize,
                *rel_threshold,
                *prefix_length,
                summary.as_ref(),
                *deplete,
                *rename,
                *threads,
                *compression_level,
                *debug,
                *quiet,
            )
            .context("Failed to run filter command")?;
        }
    }

    Ok(())
}
