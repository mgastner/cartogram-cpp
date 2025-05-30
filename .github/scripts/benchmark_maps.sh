#!/usr/bin/env bash

# Usage: benchmark_maps.sh <PR_BIN> <BASE_BIN> <OUT_JSON>

set -euo pipefail

PR_BIN=$1   # cartogram built from the pull-request
BASE_BIN=$2 # cartogram built from origin/main
OUT_JSON=$3 # file to write results to

MAP_ROOT="sample_data"
WARMUP=1   # warm-up runs
MIN_RUNS=2 # at least this many
MAX_RUNS=5 # stop early when CI width < 5% or after MAX_RUNS

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

results="$tmp/results.jsonl"
: >"$results"

run_cmd() {
  # $1 exe  $2 geo  $3 csv  $4 dir
  local exe=$1 geo=$2 csv=$3 dir=$4
  local extra=()
  [[ $dir =~ world ]] && extra+=(--world)
  "$exe" "$geo" "$csv" "${extra[@]}"
}
export -f run_cmd

for dir in "$MAP_ROOT"/*; do
  [[ -d $dir ]] || continue

  geo=$(ls "$dir"/*.geojson 2>/dev/null | head -n1) || true
  [[ -f $geo ]] || {
    echo "WARNING: $dir has no *.geojson — skipped"
    continue
  }

  for csv in "$dir"/*.csv; do
    [[ -f $csv ]] || continue
    map_name="$(basename "$dir")/$(basename "$csv")"
    echo "Benchmarking:  $map_name"

    # hyperfine runs both binaries; --ignore-failure means the overall command
    # returns 0 even if one side fails
    hyperfine --ignore-failure \
      --warmup $WARMUP \
      --min-runs $MIN_RUNS \
      --max-runs $MAX_RUNS \
      --export-json "$tmp/hf.json" \
      --command-name main \
      "bash -lc 'run_cmd $BASE_BIN \"$geo\" \"$csv\" \"$dir\"'" \
      --command-name pr \
      "bash -lc 'run_cmd $PR_BIN   \"$geo\" \"$csv\" \"$dir\"'"

    jq -n \
      --slurpfile r "$tmp/hf.json" \
      --arg map "$map_name" '
        ($r[0].results // []) as $all
        | {
            map: $map,
            base: ( ($all[]? | select(.command_name=="main" and .exit_code==0)
                      | {mean, stddev, runs}) // null ),
            pr:   ( ($all[]? | select(.command_name=="pr"   and .exit_code==0)
                      | {mean, stddev, runs}) // null )
          }' >>"$results"
  done
done

jq -s '.' "$results" >"$OUT_JSON"
echo "✅  Benchmarks written to $OUT_JSON"
