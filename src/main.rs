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
    /// Build and manipulate minimizer indexes
    Index {
        #[command(subcommand)]
        command: IndexCommands,
    },
    /// Filter fastx sequences based on presence/absence of minimizer matches
    Filter {
        /// Path to minimizer index file
        minimizers: PathBuf,

        /// Path to fastx file (or - for stdin)
        #[arg(default_value = "-")]
        input1: String,

        /// Optional path to second paired fastx file (or - for interleaved stdin)
        input2: Option<String>,

        /// Optional path to output fastx file (or - for stdout; detects .gz and .zst)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Number of minimizer matches required per query sequence (pair)
        #[arg(short = 'm', long = "matches", default_value_t = 2)]
        min_matches: usize,

        /// Consider only the first N nucleotides per sequence (0 = entire sequence)
        #[arg(short = 'n', long = "nucleotides", default_value_t = 0)]
        prefix_length: usize,

        /// Invert filtering (keep sequences WITH matches rather than those WITHOUT)
        #[arg(short = 'i', long = "invert", default_value_t = false)]
        invert: bool,

        /// Replace sequence headers with incrementing numbers
        #[arg(short = 'r', long = "rename", default_value_t = false)]
        rename: bool,

        /// Path to JSON output log file
        #[arg(short = 'l', long = "log")]
        log: Option<PathBuf>,
    },
}

#[derive(Subcommand)]
enum IndexCommands {
    /// Build index of minimizers contained within a fastx file
    Build {
        /// Path to input fastx file (supports .gz compression)
        input: PathBuf,

        /// K-mer length used for indexing
        #[arg(short = 'k', default_value_t = DEFAULT_KMER_LENGTH)]
        kmer_length: usize,

        /// Minimizer window length used for indexing
        #[arg(short = 'w', default_value_t = DEFAULT_WINDOW_SIZE)]
        window_length: usize,

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
                window_length,
                output,
            } => {
                // Convert output string to Option<PathBuf>
                let output_path = if output == "-" {
                    None
                } else {
                    Some(PathBuf::from(output))
                };

                build_index(input, *kmer_length, *window_length, output_path)
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
                output,
            } => {
                diff_index(first, second, output.as_ref())
                    .context("Failed to run index diff command")?;
            }
        },
        Commands::Filter {
            minimizers,
            input1,
            input2,
            output,
            min_matches,
            prefix_length,
            log,
            invert,
            rename,
        } => {
            run_filter(
                minimizers,
                input1,
                input2.as_deref(),
                output,
                *min_matches,
                *prefix_length,
                log.as_ref(),
                *invert,
                *rename,
            )
            .context("Failed to run filter command")?;
        }
    }

    Ok(())
}
