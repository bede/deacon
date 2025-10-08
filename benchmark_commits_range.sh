#!/bin/bash
set -e

# Benchmark all commits in a specified range
# For each commit: build index, filter single reads (2 trials), filter paired reads (2 trials)

# Configuration - Define commit range (older commit first, newer commit second)
START_COMMIT="c25fa6c"  # Older commit (exclusive)
END_COMMIT="5178d83"    # Newer commit (exclusive)

BENCH_DIR="bench"
DATA_DIR="data"
INDEX_FILE="idx.idx"

# Create bench directory if it doesn't exist
mkdir -p "$BENCH_DIR"

# Get all commits in the range (START_COMMIT..END_COMMIT)
# This gets commits from START_COMMIT (exclusive) to END_COMMIT (exclusive), in reverse chronological order
# We use ^END_COMMIT to exclude the end commit
COMMITS=$(git log --pretty=format:"%h" ${START_COMMIT}..${END_COMMIT} | grep -v "^${END_COMMIT}")

echo "Starting benchmark of commits from $START_COMMIT (exclusive) to $END_COMMIT (exclusive)..."
echo "Results will be stored in $BENCH_DIR/"

# Convert commits to array
COMMIT_ARRAY=($COMMITS)
TOTAL=${#COMMIT_ARRAY[@]}

echo "Total commits to process: $TOTAL"

CURRENT=0

for COMMIT in "${COMMIT_ARRAY[@]}"; do
    CURRENT=$((CURRENT + 1))
    echo ""
    echo "=========================================="
    echo "Processing commit $CURRENT/$TOTAL: $COMMIT"
    echo "=========================================="

    # Checkout the commit
    echo "Checking out commit $COMMIT..."
    git checkout "$COMMIT" 2>/dev/null || {
        echo "Failed to checkout $COMMIT, skipping..."
        continue
    }

    # Build the project
    echo "Building project..."
    cargo build --release > /dev/null 2>&1 || {
        echo "Build failed for $COMMIT, skipping..."
        continue
    }

    # Build index
    echo "Building index..."
    ./target/release/deacon index build "$DATA_DIR/chm13v2/chm13v2.fa.gz" > "$INDEX_FILE" 2>/dev/null || {
        echo "Index build failed for $COMMIT, skipping..."
        continue
    }

    # Single reads - Trial 1
    if [ -f "$BENCH_DIR/chm13-lr-${COMMIT}-1.json" ]; then
        echo "Single reads trial 1 already exists, skipping..."
    else
        echo "Running single reads trial 1..."
        /usr/bin/time -l ./target/release/deacon filter -d -a 1 -r 0.0 "$INDEX_FILE" \
            "$DATA_DIR/rsviruses17900/rsviruses17900.fasta" \
            -s "$BENCH_DIR/chm13-lr-${COMMIT}-1.json" > /dev/null 2>&1 || {
            echo "Single reads trial 1 failed for $COMMIT"
        }
    fi

    # Single reads - Trial 2
    if [ -f "$BENCH_DIR/chm13-lr-${COMMIT}-2.json" ]; then
        echo "Single reads trial 2 already exists, skipping..."
    else
        echo "Running single reads trial 2..."
        /usr/bin/time -l ./target/release/deacon filter -d -a 1 -r 0.0 "$INDEX_FILE" \
            "$DATA_DIR/rsviruses17900/rsviruses17900.fasta" \
            -s "$BENCH_DIR/chm13-lr-${COMMIT}-2.json" > /dev/null 2>&1 || {
            echo "Single reads trial 2 failed for $COMMIT"
        }
    fi

    # Paired reads - Trial 1
    if [ -f "$BENCH_DIR/chm13-sr-${COMMIT}-1.json" ]; then
        echo "Paired reads trial 1 already exists, skipping..."
    else
        echo "Running paired reads trial 1..."
        /usr/bin/time -l ./target/release/deacon filter -d -a 1 -r 0.0 "$INDEX_FILE" \
            "$DATA_DIR/rsviruses17900/rsviruses17900.r1.fastq.zst" \
            "$DATA_DIR/rsviruses17900/rsviruses17900.r2.fastq.zst" \
            -s "$BENCH_DIR/chm13-sr-${COMMIT}-1.json" > /dev/null 2>&1 || {
            echo "Paired reads trial 1 failed for $COMMIT"
        }
    fi

    # Paired reads - Trial 2
    if [ -f "$BENCH_DIR/chm13-sr-${COMMIT}-2.json" ]; then
        echo "Paired reads trial 2 already exists, skipping..."
    else
        echo "Running paired reads trial 2..."
        /usr/bin/time -l ./target/release/deacon filter -d -a 1 -r 0.0 "$INDEX_FILE" \
            "$DATA_DIR/rsviruses17900/rsviruses17900.r1.fastq.zst" \
            "$DATA_DIR/rsviruses17900/rsviruses17900.r2.fastq.zst" \
            -s "$BENCH_DIR/chm13-sr-${COMMIT}-2.json" > /dev/null 2>&1 || {
            echo "Paired reads trial 2 failed for $COMMIT"
        }
    fi

    echo "Completed commit $COMMIT"
done

# Return to original branch
echo ""
echo "=========================================="
echo "Benchmark complete!"
echo "Returning to main branch..."
git checkout main

echo "Results stored in $BENCH_DIR/"
echo "Total commits processed: $TOTAL"
