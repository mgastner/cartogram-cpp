#!/usr/bin/env bash

# Usage: benchmark_maps.sh <PR_BIN> <BASE_BIN> <OUT_JSON>

set -euo pipefail

PR_BIN=$1
BASE_BIN=$2
OUT_JSON=$3

MAP_ROOT="sample_data"
# WARMUP=1
# MIN_RUNS=6
# MAX_RUNS=15

WARMUP=0
MIN_RUNS=1
MAX_RUNS=1


tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

results="$tmp/results.jsonl"
: >"$results"

PY_PROCESSOR=$(
  cat <<'EOF'
import json
import sys
import os

def process_benchmark(map_name):
    try:
        with open('hf.json') as f:
            hf = json.load(f)
    except Exception:
        return {"map": map_name, "base": None, "pr": None}
    
    base = None
    pr = None
    
    for result in hf.get('results', []):
        # Skip failed runs
        if not all(code == 0 for code in result.get('exit_codes', [])):
            continue
            
        if result.get('command') == 'main':
            base = {
                'mean': result['mean'],
                'stddev': result['stddev'],
                'runs': len(result.get('times', []))
            }
        elif result.get('command') == 'pr':
            pr = {
                'mean': result['mean'],
                'stddev': result['stddev'],
                'runs': len(result.get('times', []))
            }
    
    return {
        'map': map_name,
        'base': base,
        'pr': pr
    }

if __name__ == '__main__':
    map_name = sys.argv[1]
    result = process_benchmark(map_name)
    print(json.dumps(result))
EOF
)

for dir in "$MAP_ROOT"/*; do
  [[ -d $dir ]] || continue

  geo=$(ls "$dir"/*.geojson 2>/dev/null | head -n1) || true
  [[ -f $geo ]] || {
    echo "WARNING: $dir has no *.geojson - skipped"
    continue
  }

  for csv in "$dir"/*.csv; do
    [[ -f $csv ]] || continue
    map_name="$(basename "$dir")/$(basename "$csv")"
    echo "Benchmarking:  $map_name"

    hyperfine --ignore-failure \
      --warmup "$WARMUP" \
      --min-runs "$MIN_RUNS" \
      --max-runs "$MAX_RUNS" \
      --export-json "hf.json" \
      --command-name main \
      "bash -c '[[ \$4 =~ world ]] && set -- \"\$1\" \"\$2\" \"\$3\" --world || set -- \"\$1\" \"\$2\" \"\$3\"; exec \"\$@\"' _ \
         \"$BASE_BIN\" \"$geo\" \"$csv\" \"$dir\"" \
      --command-name pr \
      "bash -c '[[ \$4 =~ world ]] && set -- \"\$1\" \"\$2\" \"\$3\" --world || set -- \"\$1\" \"\$2\" \"\$3\"; exec \"\$@\"' _ \
         \"$PR_BIN\"  \"$geo\" \"$csv\" \"$dir\""

    python3 -c "$PY_PROCESSOR" "$map_name" >>"$results"
  done
done

# Final aggregation
python3 -c "import json, sys; data = [json.loads(line) for line in sys.stdin]; print(json.dumps(data))" <"$results" >"$OUT_JSON"

echo "✅  Benchmarks written to $OUT_JSON"
