#!/usr/bin/env python3
"""Run the conformality experiments of the ICC 2027 paper and tabulate them.

Each run computes the cartogram T0 of one map with the flow-based method,
without polygon simplification (-P 100000000), and then the correction A
with the production settings: alpha = 1/2, ocean weight 1, L-BFGS memory 8,
Stormer-Verlet with substep factor 1, and evaluation quadrature 4. The runs
differ in the map and in the number of modes per side.

So that every number can be traced to a commit, the script
  1. refuses to run if the sources of the repository have uncommitted or
     untracked changes (unless --allow-dirty),
  2. builds the cartogram target in BUILD_DIR, so that the binary matches
     the commit,
  3. writes provenance.json into each run directory: the commit, the SHA-256
     of the binary, the command, the number of OpenMP threads, the host,
     the exit code, and the wall time,
  4. writes OUT_DIR/tables.md: the provenance of the runs, followed by the
     tables of summarize_conformality.py.
A run whose directory already has a provenance.json with the same commit
and binary, and exit code 0, is not repeated (use --force to repeat it).

On macOS, caffeinate keeps the machine awake while the script runs.

Only the Python standard library is used.

Usage: paper_runs.py BUILD_DIR OUT_DIR [--threads N] [--only NAME[,NAME...]]
                     [--allow-dirty] [--force]
"""

import argparse
import datetime
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SAMPLE_DATA = REPO / "sample_data"

MAPS = {
    "usa": (
        "conterminous_usa_by_state_since_1959/"
        "conterminous_usa_by_state_since_1959.geojson",
        "conterminous_usa_by_state_since_1959/"
        "conterminous_usa_population_by_state_2020.csv",
    ),
    "china": (
        "china_by_province_with_chinese_taipei_since_1998/"
        "china_by_province_with_chinese_taipei_since_1998.geojson",
        "china_by_province_with_chinese_taipei_since_1998/"
        "china_population_by_province_2020.csv",
    ),
    "sui": (
        "switzerland_by_canton_since_1999/"
        "switzerland_by_canton_since_1999.geojson",
        "switzerland_by_canton_since_1999/"
        "switzerland_population_by_canton_2016.csv",
    ),
    "world": (
        "world_by_country_since_2022/world_by_country_since_2022.geojson",
        "world_by_country_since_2022/world_population_by_country_2010.csv",
    ),
}

# Run name, map, modes per side
RUNS = [
    ("usa16", "usa", 16),
    ("usa32", "usa", 32),
    ("usa64", "usa", 64),
    ("china16", "china", 16),
    ("china32", "china", 32),
    ("china64", "china", 64),
    ("sui32", "sui", 32),
    ("sui64", "sui", 64),
    ("world32", "world", 32),
    ("world64", "world", 64),
]

SETTINGS = [
    "-P", "100000000",
    "--conf_optimize",
    "--conf_alpha", "0.5",
    "--conf_lbfgs", "8",
    "--conf_integrator", "verlet",
    "--conf_substep_factor", "1",
    "--conf_eval_quad", "4",
]


def git(*args):
    return subprocess.run(
        ["git", "-C", str(REPO), *args],
        check=True, capture_output=True, text=True,
    ).stdout.strip()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def cmake_command(build_dir):
    # Build with the CMake that configured the build directory
    cache = Path(build_dir) / "CMakeCache.txt"
    m = re.search(r"^CMAKE_COMMAND:INTERNAL=(.+)$", cache.read_text(), re.M)
    if not m:
        sys.exit(f"{cache}: no CMAKE_COMMAND")
    return m.group(1)


def main():
    p = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("build_dir", type=Path)
    p.add_argument("out_dir", type=Path)
    p.add_argument("--threads", type=int, default=os.cpu_count())
    p.add_argument("--only", help="comma-separated run names")
    p.add_argument("--allow-dirty", action="store_true")
    p.add_argument("--force", action="store_true")
    args = p.parse_args()

    runs = RUNS
    if args.only:
        names = args.only.split(",")
        unknown = set(names) - {r[0] for r in RUNS}
        if unknown:
            sys.exit(f"unknown runs: {', '.join(sorted(unknown))}")
        runs = [r for r in RUNS if r[0] in names]

    # CMake compiles every source file under src/, tracked or not
    dirty = git(
        "status", "--porcelain", "--", "src", "include", "CMakeLists.txt",
        "conanfile.txt",
    )
    if dirty and not args.allow_dirty:
        sys.exit(f"uncommitted changes in the sources:\n{dirty}")
    commit = git("rev-parse", "HEAD")
    if dirty:
        commit += "-dirty"

    subprocess.run(
        [cmake_command(args.build_dir), "--build", str(args.build_dir),
         "--target", "cartogram"],
        check=True,
    )
    binary = (args.build_dir / "cartogram").resolve()
    binary_sha = sha256(binary)

    if sys.platform == "darwin" and shutil.which("caffeinate"):
        subprocess.Popen(["caffeinate", "-i", "-w", str(os.getpid())])

    args.out_dir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads))
    for name, map_name, modes in runs:
        run_dir = args.out_dir / name
        prov_file = run_dir / "provenance.json"
        if prov_file.exists() and not args.force:
            old = json.loads(prov_file.read_text())
            if (old["commit"], old["binary_sha256"], old["exit_code"]) == (
                commit, binary_sha, 0
            ):
                print(f"{name}: done before, skipped", flush=True)
                continue
        run_dir.mkdir(exist_ok=True)
        geo, csv = (str(SAMPLE_DATA / f) for f in MAPS[map_name])
        cmd = [str(binary), geo, csv, *SETTINGS, "--conf_modes", str(modes)]
        if map_name == "world":
            cmd.append("--world")
        start = datetime.datetime.now().astimezone()
        t0 = time.perf_counter()
        with open(run_dir / "out.log", "w") as out, \
                open(run_dir / "err.log", "w") as err:
            code = subprocess.run(
                cmd, cwd=run_dir, stdout=out, stderr=err, env=env
            ).returncode
        wall = time.perf_counter() - t0
        prov = {
            "commit": commit,
            "binary": str(binary),
            "binary_sha256": binary_sha,
            "command": cmd,
            "omp_num_threads": args.threads,
            "host": platform.node(),
            "platform": platform.platform(),
            "start": start.isoformat(timespec="seconds"),
            "wall_seconds": round(wall, 2),
            "exit_code": code,
        }
        prov_file.write_text(json.dumps(prov, indent=2) + "\n")
        print(f"{name}: exit {code} after {wall:.0f} s", flush=True)

    # Tables of all runs in OUT_DIR that finished, in the order of RUNS
    done = []
    for name, _, _ in RUNS:
        prov_file = args.out_dir / name / "provenance.json"
        if prov_file.exists():
            prov = json.loads(prov_file.read_text())
            if prov["exit_code"] == 0:
                done.append((name, prov))
    if not done:
        sys.exit("no finished runs")
    lines = [
        "# Conformality runs of the ICC 2027 paper",
        "",
        "| run | commit | binary SHA-256 | threads | start | wall time (s) |",
        "|---|---|---|---|---|---|",
    ]
    for name, prov in done:
        lines.append(
            f"| {name} | {prov['commit'][:12]} | {prov['binary_sha256'][:12]} "
            f"| {prov['omp_num_threads']} | {prov['start']} "
            f"| {prov['wall_seconds']:.0f} |"
        )
    summary = subprocess.run(
        [sys.executable, str(Path(__file__).with_name(
            "summarize_conformality.py"))]
        + [str(args.out_dir / name) for name, _ in done],
        check=True, capture_output=True, text=True,
    ).stdout
    (args.out_dir / "tables.md").write_text(
        "\n".join(lines) + "\n\n" + summary
    )
    print(f"wrote {args.out_dir / 'tables.md'}")


if __name__ == "__main__":
    main()
