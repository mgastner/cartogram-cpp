#!/usr/bin/env python3
import argparse
import csv
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time

sys.stdout.reconfigure(line_buffering=True)


def parse_args():
    p = argparse.ArgumentParser(description="Cartogram per-map fuzz-test harness")
    p.add_argument("data_dir", help="Path to data directory")
    p.add_argument("cartogram_bin", help="Path to cartogram executable")
    p.add_argument(
        "--max-rows",
        type=int,
        default=100,
        help="Skip maps that have more regions than this",
    )
    p.add_argument(
        "--runs", type=int, default=10, help="Number of fuzz iterations per map"
    )
    p.add_argument(
        "--out-dir", default="fuzzer_csvs", help="Base directory to save failing CSVs"
    )
    p.add_argument(
        "--map",
        dest="single_map",
        required=True,
        help="Name of the map (subfolder under sample_data) to fuzz",
    )
    return p.parse_args()


def generate_value():
    if random.random() < 0.5:
        return random.random()
    else:
        exp = random.uniform(3, 6)
        return 10**exp


def main():
    args = parse_args()
    sample_root = args.data_dir
    carto_bin = args.cartogram_bin
    max_rows = args.max_rows
    runs = args.runs

    # Directory to save failing CSVs
    fail_dir = os.path.abspath(os.path.join(args.out_dir, args.single_map))

    # Clear only this map's past run failures, so other maps remain intact
    if os.path.isdir(fail_dir):
        shutil.rmtree(fail_dir)
    os.makedirs(fail_dir, exist_ok=True)

    # Locate the map folder
    dirp = os.path.join(sample_root, args.single_map)
    if not os.path.isdir(dirp):
        print(f"[ERROR] map folder not found: {dirp}", file=sys.stderr)
        return 1

    # Find the geojson
    geojsons = [f for f in os.listdir(dirp) if f.lower().endswith(".geojson")]
    if not geojsons:
        print(f"[SKIP] no .geojson in {dirp}", file=sys.stderr)
        return 77

    geo_path = os.path.join(dirp, geojsons[0])

    # Find all CSVs
    csvs = sorted(f for f in os.listdir(dirp) if f.lower().endswith(".csv"))

    if not csvs:
        print(f"[SKIP] no .csv in {dirp}", file=sys.stderr)
        return 77

    # Pick any csv
    csv_name = csvs[0]
    csv_path = os.path.join(dirp, csv_name)

    with open(csv_path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = list(csv.reader(fh))

    if len(reader) - 1 > max_rows:
        print(f"[SKIP] {csv_name}: {len(reader)} rows > {max_rows}")
        return 77

    header = reader[0]
    data_rows = reader[1:]
    print(f"\n=== Map: {args.single_map} | CSV: {csv_name} ({len(data_rows)} rows) ===")

    n_fail = 0
    for i in range(1, runs + 1):
        # Prepare the fuzzed CSV
        fuzzed = [header]
        for row in data_rows:
            val0 = row[0]
            val1 = f"{generate_value():.6f}"
            remainder = row[2:] if len(row) > 2 else []
            fuzzed.append([val0, val1] + remainder)

        fd, tmp_csv = tempfile.mkstemp(
            prefix=f"fuzz_{args.single_map}_", suffix="_run{i}.csv"
        )
        os.close(fd)
        with open(tmp_csv, "w", newline="") as outfh:
            csv.writer(outfh).writerows(fuzzed)

        # Run
        cmd = [carto_bin, geo_path, tmp_csv]

        # Make sure the cartogram runs fast
        cmd.extend(["--max_allowed_autoscale_grid_length", "512"])

        if "world" in args.single_map.lower():
            cmd.append("--world")

        start = time.perf_counter()
        proc = subprocess.run(cmd, capture_output=True, text=True)
        dur = time.perf_counter() - start

        ok = proc.returncode == 0
        tag = "[ OK ]" if ok else "[FAIL]"
        print(f"    {tag} run {i:2d}/{runs} in {dur:.2f}s")

        fail_name = f"{args.single_map}_run{i}.csv"
        fail_path = os.path.join(fail_dir, fail_name)
        stderr_path = os.path.splitext(fail_path)[0] + ".stderr"

        if ok:
            if os.path.exists(fail_path):
                os.remove(fail_path)
            if os.path.exists(stderr_path):
                os.remove(stderr_path)
        else:
            n_fail += 1
            shutil.copy(tmp_csv, fail_path)

            with open(stderr_path, "w") as errfh:
                errfh.write(proc.stderr)

        os.remove(tmp_csv)

    print("\nFuzzer summary:", "FAILURES detected." if n_fail else "All passed.")
    print()
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
