#!/usr/bin/env python3
"""Write the tables of the ICC 2027 paper as LaTeX, from the outputs of paper_runs.py.

Usage: paper_tables_tex.py OUT_DIR [--output FILE] [--timing DIR ...]

OUT_DIR is a directory written by paper_runs.py: one subdirectory per run,
each with provenance.json, err.log, and the conformality JSON files. The
script writes three tables (booktabs, as loaded by the paper):

  T1  T0 against T at the finest mode count of each map, for land, ocean,
      and the whole rectangle (evaluation cloud).
  T2  verification checks of every run.
  T3  cost and mode dependence of every run.

No number is typed by hand: every value comes from the JSON files, the time
report in err.log, or provenance.json, and the table captions name OUT_DIR's
commit and thread count.
"""

import argparse
import glob
import json
import math
import os
import re
import sys

MAP_NAMES = {"usa": "USA", "china": "China", "sui": "Switzerland", "world": "World"}
MAP_ORDER = ["usa", "china", "sui", "world"]


def load_run(run_dir):
    """Returns a dict with the JSON outputs and timings of one run."""
    prov = json.load(open(os.path.join(run_dir, "provenance.json")))
    opt_files = glob.glob(os.path.join(run_dir, "*_conformality_opt.json"))
    t0_files = glob.glob(os.path.join(run_dir, "*_conformality_t0.json"))
    if len(opt_files) != 1 or len(t0_files) != 1:
        sys.exit(f"{run_dir}: expected one *_opt.json and one *_t0.json")
    opt = json.load(open(opt_files[0]))
    t0 = json.load(open(t0_files[0]))
    err = open(os.path.join(run_dir, "err.log"), errors="replace").read()
    m = re.search(r"^Integration: (\d+) ms$", err, re.M)
    if not m:
        sys.exit(f"{run_dir}: no 'Integration' line in the time report")
    name = os.path.basename(run_dir.rstrip("/"))
    key = re.match(r"([a-z]+)(\d+)$", name)
    if not key:
        sys.exit(f"{run_dir}: run name must be <map><modes>, e.g. china64")
    return {
        "name": name,
        "dir": run_dir,
        "map": key.group(1),
        "modes_long": int(key.group(2)),
        "prov": prov,
        "opt": opt,
        "t0": t0,
        "t0_seconds": int(m.group(1)) / 1000.0,
    }


def sig(x, digits=3):
    """Formats x with the given number of significant digits."""
    if x == 0 or not math.isfinite(x):
        return f"{x:g}"
    e = math.floor(math.log10(abs(x)))
    decimals = max(0, digits - 1 - e)
    return f"{x:.{decimals}f}"


def sci(x):
    """Formats a small number as m \\times 10^{e} for LaTeX."""
    if x == 0:
        return "0"
    e = math.floor(math.log10(abs(x)))
    m = x / 10**e
    return f"${m:.1f}\\times10^{{{e}}}$"


def change(before, after):
    return f"{100.0 * (after / before - 1.0):+.1f}\\%".replace("-", "$-$")


def arrow(before, after):
    return f"{sig(before)} $\\to$ {sig(after)}"


def modes_label(run):
    m = run["opt"]["modes"]
    return f"${m[0]}\\times{m[1]}$"


def finest_runs(runs):
    best = {}
    for r in runs:
        if r["map"] not in best or r["modes_long"] > best[r["map"]]["modes_long"]:
            best[r["map"]] = r
    return [best[k] for k in MAP_ORDER if k in best]


def table_t1(runs, provenance_note):
    rows = []
    for r in finest_runs(runs):
        b, a = r["opt"]["before"], r["opt"]["after"]
        for subset, label in (("land", "land"), ("ocean", "ocean"), ("rectangle", "all")):
            suffix = "" if subset == "rectangle" else f"_{subset}"
            e0b = b["evaluation_cloud"]["E0" + suffix]
            e0a = a["evaluation_cloud"]["E0" + suffix]
            ehb = b["evaluation_cloud"]["E_half" + suffix]
            eha = a["evaluation_cloud"]["E_half" + suffix]
            db, da = b[f"distortion_{subset}"], a[f"distortion_{subset}"]
            first = label == "land"
            rows.append(
                (f"{MAP_NAMES[r['map']]} ({modes_label(r)})" if first else "")
                + f" & {label}"
                + f" & {arrow(e0b, e0a)} & {change(e0b, e0a)}"
                + f" & {arrow(ehb, eha)} & {change(ehb, eha)}"
                + f" & {arrow(db['e_mean'], da['e_mean'])}"
                + f" & {arrow(db['e_p90'], da['e_p90'])}"
                + f" & {arrow(db['e_p99'], da['e_p99'])} \\\\"
            )
        rows.append("\\addlinespace")
    rows = rows[:-1]
    return "\n".join(
        [
            "\\begin{table*}",
            "\\centering",
            "\\caption{Distortion of the cartogram $T_0$ and of the corrected map $T$ at the finest mode count of each map, "
            "on the evaluation cloud. $E_0$ and $E_{1/2}$ are means of the excess MIPS energy; $e = \\ln(a/b)$, with mean, "
            f"90th and 99th percentile. {provenance_note}}}",
            "\\label{tab:effect}",
            "\\footnotesize",
            "\\setlength{\\tabcolsep}{3.5pt}",
            "\\begin{tabular}{llrrrrrrr}",
            "\\toprule",
            "Map & Subset & \\multicolumn{2}{c}{$E_0$} & \\multicolumn{2}{c}{$E_{1/2}$} & mean $e$ & $e$, p90 & $e$, p99 \\\\",
            "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
            " & & $T_0 \\to T$ & change & $T_0 \\to T$ & change & $T_0 \\to T$ & $T_0 \\to T$ & $T_0 \\to T$ \\\\",
            "\\midrule",
            *rows,
            "\\bottomrule",
            "\\end{tabular}",
            "\\end{table*}",
        ]
    )


def mesh_check_file(run_dir):
    """Ratio and error of the two finest cell refinements in mesh_check_*.txt,
    written by conformality_check, or None if the run has no such file."""
    files = sorted(glob.glob(os.path.join(run_dir, "mesh_check_*.txt")))
    if not files:
        return None
    rows = []
    for line in open(files[-1]):
        f = line.split()
        if len(f) >= 4 and f[0].isdigit() and f[1] == "1":
            rows.append((int(f[0]), float(f[3])))
    if len(rows) < 2:
        return None
    rows.sort()
    return rows[-2][1] / rows[-1][1], rows[-1][0]


def finest_ratio(entries, key_size, key_err):
    """Convergence ratio of the two finest refinements, and the finest error."""
    entries = sorted(entries, key=lambda e: e[key_size])
    if len(entries) < 2:
        return None, None
    return entries[-2][key_err] / entries[-1][key_err], entries[-1][key_err]


def replay_cell(c):
    # Largest difference between a replay of the step log and the optimized
    # nodes (positions and Jacobians); 0 means bit for bit
    if not c["replay_converged"]:
        return "failed"
    m = c.get("replay_mismatch_optimization_cloud")
    return "--" if m is None else ("0" if m == 0 else sci(m))


def table_t2(runs, provenance_note):
    rows = []
    for r in runs:
        c = r["opt"]["checks"]
        bd = c["boundary"]
        boundary = max(bd["max_corner_motion"], bd["max_offset_from_side"], bd["max_overshoot_along_side"])
        fd = c["finite_differences"]
        from_file = mesh_check_file(r["dir"])
        if from_file:
            ratio, finest = from_file
        else:
            mesh_cells = [e for e in c.get("mesh", []) if e["samples_per_edge"] == 1]
            ratio, finest = finest_ratio(mesh_cells, "cells_long", "max_rel_area_error")
            finest = max((e["cells_long"] for e in mesh_cells), default=None)
        areas = c.get("region_areas", [])
        non_simple = sum(e["non_simple_rings"] for e in areas)
        area_err = min((e["max_rel_area_error"] for e in areas), default=float("nan"))
        n_nonpos = r["opt"]["after"]["evaluation_cloud"]["n_non_positive_det"]
        rows.append(
            f"{MAP_NAMES[r['map']]} ({modes_label(r)}) & {sci(c['max_det_drift_evaluation_cloud'])}"
            f" & {n_nonpos} & {sci(boundary) if boundary else '0'}"
            f" & {sci(fd['richardson_max_rel_error'])}"
            f" & {replay_cell(c)}"
            f" & {sig(ratio, 3) if ratio else '--'} ({finest if finest else '--'}) & {sci(area_err)} & {non_simple} \\\\"
        )
    return "\n".join(
        [
            "\\begin{table*}",
            "\\centering",
            "\\caption{Verification of every run: largest drift of $\\det DT$ from $\\det DT_0$ on the evaluation cloud, "
            "number of nodes with $\\det DT \\le 0$, largest motion of the boundary, largest relative error of the "
            "propagated Jacobian $DT$ against finite differences of tracer positions (Richardson extrapolation, "
            "1001 nodes), largest difference between a replay of the step log and the optimized nodes "
            "(0: bit for bit), ratio of the "
            "largest cell-area errors at the two finest mesh refinements, with the number of cells along the longer side "
            "of the finest one in parentheses (4 for second-order convergence), largest "
            "relative area error of a corrected region at the finest polygon resolution, and number of rings that are "
            f"not simple. {provenance_note}}}",
            "\\label{tab:verification}",
            "\\footnotesize",
            "\\setlength{\\tabcolsep}{3.5pt}",
            "\\begin{tabular}{lrrrrrrrr}",
            "\\toprule",
            "Run & det drift & $\\det \\le 0$ & boundary & $DT$ check & replay & mesh ratio & area error & non-simple \\\\",
            "\\midrule",
            *rows,
            "\\bottomrule",
            "\\end{tabular}",
            "\\end{table*}",
        ]
    )


def median(values):
    v = sorted(values)
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def timings(run, timing_sets):
    """Median wall times of a run over the timing sets (or, without timing
    sets, of the run itself), and the range of its total time."""
    samples = [t[run["name"]] for t in timing_sets if run["name"] in t] or [run]
    t0 = median([s["t0_seconds"] for s in samples])
    opt = median([s["opt"]["seconds"]["optimize"] for s in samples])
    replay = median([s["opt"]["seconds"]["replay"] for s in samples])
    totals = [s["prov"]["wall_seconds"] for s in samples]
    return t0, opt, replay, median(totals), min(totals), max(totals), len(samples)


def table_t3(runs, provenance_note, timing_sets=()):
    rows = []
    previous_map = None
    n_samples = set()
    threads = set()
    for r in runs:
        if previous_map is not None and r["map"] != previous_map:
            rows.append("\\addlinespace")
        previous_map = r["map"]
        o = r["opt"]
        ev_b, ev_a = o["before"]["evaluation_cloud"], o["after"]["evaluation_cloud"]
        lat = o["lattice"]
        t0, opt, replay, total, lo, hi, n = timings(r, timing_sets)
        n_samples.add(n)
        timed = [t[r["name"]] for t in timing_sets if r["name"] in t] or [r]
        threads.update(t["prov"]["omp_num_threads"] for t in timed)
        rows.append(
            f"{MAP_NAMES[r['map']]} & ${lat[0]}\\times{lat[1]}$ & {modes_label(r)} & {o['n_accepted_steps']}"
            f" & {t0:.1f} & {opt:.0f} & {replay:.0f}"
            f" & {total:.0f} ({lo:.0f}--{hi:.0f})"
            f" & {change(ev_b['E_half'], ev_a['E_half'])} & {change(ev_b['E0_land'], ev_a['E0_land'])} \\\\"
        )
    return "\n".join(
        [
            "\\begin{table*}",
            "\\centering",
            "\\caption{Cost and mode dependence: lattice of $T_0$, number of modes, accepted steps, wall times in seconds "
            "(integration of $T_0$, optimization, application of the correction to the evaluation cloud and "
            "the polygons together with the area and mesh checks, whole run with its range in "
            "parentheses), and the change of $E_{1/2}$ "
            f"over the rectangle and of $E_0$ on land. Times are medians of {'/'.join(map(str, sorted(n_samples)))} runs "
            f"with {'/'.join(map(str, sorted(threads)))} threads; the other columns come from the runs named below. "
            f"{provenance_note}}}",
            "\\label{tab:cost}",
            "\\footnotesize",
            "\\setlength{\\tabcolsep}{3.5pt}",
            "\\begin{tabular}{lllrrrrrrr}",
            "\\toprule",
            "Map & Lattice & Modes & Steps & $T_0$ & optimize & replay & total & $E_{1/2}$ & $E_0$, land \\\\",
            "\\midrule",
            *rows,
            "\\bottomrule",
            "\\end{tabular}",
            "\\end{table*}",
        ]
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out_dir")
    parser.add_argument("--output", help="LaTeX file to write (default: stdout)")
    parser.add_argument(
        "--timing", nargs="*", default=[],
        help="more output directories of paper_runs.py whose wall times enter the medians of T3; "
        "their outputs must equal those of OUT_DIR",
    )
    args = parser.parse_args()

    run_dirs = sorted(d for d in glob.glob(os.path.join(args.out_dir, "*")) if os.path.isfile(os.path.join(d, "provenance.json")))
    if not run_dirs:
        sys.exit(f"{args.out_dir}: no run directories with provenance.json")
    runs = [load_run(d) for d in run_dirs]
    runs.sort(key=lambda r: (MAP_ORDER.index(r["map"]), r["modes_long"]))

    commits = {r["prov"]["commit"][:12] for r in runs}
    threads = {r["prov"]["omp_num_threads"] for r in runs}
    failed = [r["name"] for r in runs if r["prov"]["exit_code"] != 0]
    if failed:
        sys.exit(f"runs with a nonzero exit code: {', '.join(failed)}")
    note = (
        f"Commit {', '.join(sorted(commits))}; {', '.join(str(t) for t in sorted(threads))} threads; "
        f"$\\alpha = {runs[0]['opt']['settings']['alpha']:g}$."
    )

    timing_sets = []
    for d in args.timing:
        tset = {}
        for rd in sorted(glob.glob(os.path.join(d, "*"))):
            if os.path.isfile(os.path.join(rd, "provenance.json")):
                t = load_run(rd)
                tset[t["name"]] = t
        for r in runs:
            t = tset.get(r["name"])
            if t and t["opt"]["after"]["evaluation_cloud"] != r["opt"]["after"]["evaluation_cloud"]:
                sys.exit(f"{d}/{r['name']}: the result differs from {args.out_dir}; not a timing repeat")
        timing_sets.append(tset)

    text = "\n\n".join(
        [
            f"% Generated by research/scripts/paper_tables_tex.py from {os.path.abspath(args.out_dir)}",
            table_t1(runs, note),
            table_t2(runs, note),
            table_t3(runs, note, timing_sets),
        ]
    ) + "\n"
    if args.output:
        open(args.output, "w").write(text)
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
