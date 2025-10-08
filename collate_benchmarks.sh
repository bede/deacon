#!/bin/bash

# Collate bp_per_second from all benchmark JSON files
# Output format: commit, type (lr/sr), trial1, trial2, max

BENCH_DIR="bench"
OUTPUT_FILE="benchmark_results.tsv"

echo "Collating benchmark results from $BENCH_DIR..."

# Write header
echo -e "commit\ttype\ttrial1_bp_per_second\ttrial2_bp_per_second\tmax_bp_per_second" > "$OUTPUT_FILE"

# Get unique commits from JSON files
UNIQUE_COMMITS=$(for json_file in "$BENCH_DIR"/*.json; do
    [ -e "$json_file" ] || continue
    filename=$(basename "$json_file")
    if [[ $filename =~ chm13-(lr|sr)-([^-]+)-([0-9]+)\.json ]]; then
        echo "${BASH_REMATCH[2]}"
    fi
done | sort -u)

# Get commit order from git history
COMMIT_ORDER=$(git log --all --pretty=format:"%h" | grep -F -f <(echo "$UNIQUE_COMMITS"))

# Process commits in reverse chronological order
for commit in $COMMIT_ORDER; do
    # Process both lr and sr types for this commit
    for type in lr sr; do
        # Check if this commit-type combination exists
        if [ ! -f "$BENCH_DIR/chm13-${type}-${commit}-1.json" ] && [ ! -f "$BENCH_DIR/chm13-${type}-${commit}-2.json" ]; then
            continue
        fi
        # Get bp_per_second from trial 1
        trial1_file="$BENCH_DIR/chm13-${type}-${commit}-1.json"
        trial2_file="$BENCH_DIR/chm13-${type}-${commit}-2.json"

        if [ -f "$trial1_file" ]; then
            trial1_bp=$(grep -o '"bp_per_second":[^,}]*' "$trial1_file" | cut -d: -f2 | tr -d ' ')
        else
            trial1_bp="NA"
        fi

        if [ -f "$trial2_file" ]; then
            trial2_bp=$(grep -o '"bp_per_second":[^,}]*' "$trial2_file" | cut -d: -f2 | tr -d ' ')
        else
            trial2_bp="NA"
        fi

        # Calculate max
        if [ "$trial1_bp" != "NA" ] && [ "$trial2_bp" != "NA" ]; then
            max_bp=$(echo "$trial1_bp $trial2_bp" | awk '{print ($1 > $2) ? $1 : $2}')
        elif [ "$trial1_bp" != "NA" ]; then
            max_bp="$trial1_bp"
        elif [ "$trial2_bp" != "NA" ]; then
            max_bp="$trial2_bp"
        else
            max_bp="NA"
        fi

        # Output row
        echo -e "${commit}\t${type}\t${trial1_bp}\t${trial2_bp}\t${max_bp}" >> "$OUTPUT_FILE"
    done
done

echo "Results written to $OUTPUT_FILE"
echo ""
echo "Summary:"
cat "$OUTPUT_FILE" | column -t -s $'\t'
