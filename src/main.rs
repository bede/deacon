mod filter;
mod index;
mod index_build;
mod index_diff;
mod index_format;
mod index_info;
mod index_union;
mod minimizers;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build and manipulate minimizer indexes
    Index {
        #[command(subcommand)]
        command: IndexCommands,
    },
    /// Filter fastx sequences based on presence/absence of minimizer matches
    Filter {
        /// Minimizer index file path
        minimizers: PathBuf,

        /// Input fastx file path (- for stdin)
        #[arg(default_value = "-")]
        input: String,

        /// Path to output fastx file (- for stdout; detects .gz and .zst extensions)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Minimum number of minimizer matches per query sequence
        #[arg(short = 'm', long = "matches", default_value_t = 1)]
        min_matches: usize,

        /// Consider only the first N nucleotides per sequence (0 = entire sequence)
        #[arg(short = 'n', long = "nucleotides", default_value_t = 0)]
        prefix_length: usize,

        /// Path to JSON report file
        #[arg(long = "report")]
        report: Option<PathBuf>,

        /// Invert filtering (keep sequences WITH matches rather than those WITHOUT)
        #[arg(short = 'i', long = "invert", default_value_t = false)]
        invert: bool,

        /// Replace sequence headers with sequential numbers (1, 2, 3...)
        #[arg(short = 'r', long = "rename", default_value_t = false)]
        rename: bool,
    },
}

#[derive(Subcommand)]
enum IndexCommands {
    /// Build index of minimizers contained within a FASTX file
    Build {
        /// Path to input FASTX file (supports .gz compression)
        input: PathBuf,

        /// K-mer length used for indexing
        #[arg(short = 'k', default_value_t = minimizers::DEFAULT_KMER_LENGTH)]
        kmer_length: usize,

        /// Window size used for indexing
        #[arg(short = 'w', default_value_t = minimizers::DEFAULT_WINDOW_SIZE)]
        window_size: usize,

        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,
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

        /// Path to second index file (to subtract from first)
        #[arg(required = true)]
        second: PathBuf,

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
            } => {
                // Convert output string to Option<PathBuf>
                let output_path = if output == "-" {
                    None
                } else {
                    Some(PathBuf::from(output))
                };

                index::build(input, *kmer_length, *window_size, output_path)
                    .context("Failed to run index build command")?;
            }
            IndexCommands::Info { index } => {
                index::info(index).context("Failed to run index info command")?;
            }
            IndexCommands::Union { inputs, output } => {
                index::union(inputs, output.as_ref())
                    .context("Failed to run index union command")?;
            }
            IndexCommands::Diff {
                first,
                second: subtract,
                output,
            } => {
                index::diff(first, subtract, output.as_ref())
                    .context("Failed to run index diff command")?;
            }
        },
        Commands::Filter {
            minimizers,
            input,
            output,
            min_matches,
            prefix_length,
            report,
            invert,
            rename,
        } => {
            filter::run(
                minimizers,
                input,
                output,
                *min_matches,
                *prefix_length,
                report.as_ref(),
                *invert,
                *rename,
            )
            .context("Failed to run filter command")?;
        }
    }

    Ok(())
}
