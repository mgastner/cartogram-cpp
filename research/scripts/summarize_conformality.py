#!/usr/bin/env python3
"""Summarize conformality optimizer runs as Markdown tables.

Each RUN_DIR must contain one <inset>_conformality_opt.json written by
cartogram with --conf_optimize. Prints these tables to stdout:

  1. effect sizes (E0 and e = ln(a/b)) on the evaluation cloud,
  2. both weightings of the energy, E0 (area on the input map) and E1
     (area on the cartogram), and which one the run minimized (alpha), for
     runs that report E1,
  3. distortion by subset (rectangle, land, ocean), for runs that report
     ocean statistics,
  4. correctness checks (det drift, check mesh, region areas, FD, boundary),
  5. quadrature agreement (optimization versus evaluation cloud) and cost.

Only the Python standard library is used.

Usage: summarize_conformality.py RUN_DIR [RUN_DIR ...]
"""

import json
import sys
from pathlib import Path


def load(run_dir):
    files = sorted(Path(run_dir).glob("*_conformality_opt.json"))
    if len(files) != 1:
        sys.exit(f"{run_dir}: expected one *_conformality_opt.json")
    return json.loads(files[0].read_text())


def change(before, after):
    return f"{before:.4g} → {after:.4g} ({100 * (after / before - 1):+.1f}%)"


def ratios(values):
    return ", ".join(f"{a / b:.1f}" for a, b in zip(values, values[1:]) if b > 0)


def table(header, rows):
    lines = ["| " + " | ".join(header) + " |"]
    lines.append("|" + "---|" * len(header))
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines)


def main():
    runs = [(Path(d).name, load(d)) for d in sys.argv[1:]]
    if not runs:
        sys.exit(__doc__)

    effect = []
    weighting = []
    subsets = []
    checks = []
    quadrature = []
    for name, r in runs:
        b, a = r["before"], r["after"]
        be, ae = b["evaluation_cloud"], a["evaluation_cloud"]
        bl, al = b["distortion_land"], a["distortion_land"]
        br, ar = b["distortion_rectangle"], a["distortion_rectangle"]
        modes = "×".join(str(m) for m in r["modes"])
        # Runs before the ocean statistics were added lack these fields
        has_ocean = "E0_ocean" in be and "distortion_ocean" in b
        effect.append(
            [
                name,
                modes,
                change(be["E0"], ae["E0"]),
                change(be["E0_land"], ae["E0_land"]),
                change(be["E0_ocean"], ae["E0_ocean"]) if has_ocean else "n/a",
                change(br["e_mean"], ar["e_mean"]),
                change(bl["e_mean"], al["e_mean"]),
                change(b["distortion_ocean"]["e_mean"], a["distortion_ocean"]["e_mean"])
                if has_ocean
                else "n/a",
                change(bl["e_p99"], al["e_p99"]),
            ]
        )
        if "E1" in be:
            weighting.append(
                [name, modes, r["settings"].get("alpha", 0)]
                + [
                    change(be[f"{e}{s}"], ae[f"{e}{s}"])
                    for s in ("", "_land", "_ocean")
                    for e in ("E0", "E1")
                ]
            )
        if has_ocean:
            for subset, energy_key in (
                ("rectangle", "E0"),
                ("land", "E0_land"),
                ("ocean", "E0_ocean"),
            ):
                bd, ad = b[f"distortion_{subset}"], a[f"distortion_{subset}"]
                subsets.append(
                    [
                        name,
                        subset,
                        f"{100 * bd['n_nodes'] / be['n_nodes']:.1f}%",
                        change(be[energy_key], ae[energy_key]),
                        change(bd["e_mean"], ad["e_mean"]),
                        change(bd["angular_mean"], ad["angular_mean"]),
                        change(bd["e_p50"], ad["e_p50"]),
                        change(bd["e_p90"], ad["e_p90"]),
                        change(bd["e_p99"], ad["e_p99"]),
                        change(bd["angular_p99"], ad["angular_p99"]),
                    ]
                )

        c = r["checks"]
        mesh = c["mesh"]
        straight = [m["max_rel_area_error"] for m in mesh if m["samples_per_edge"] == 1]
        curved = [mesh[0]["max_rel_area_error"]] + [
            m["max_rel_area_error"] for m in mesh if m["samples_per_edge"] > 1
        ]
        regions = c["region_areas"]
        fd = c["finite_differences"]
        checks.append(
            [
                name,
                f"{max(c['max_det_drift_optimization_cloud'], c['max_det_drift_evaluation_cloud']):.1e}",
                f"{straight[-1]:.1e} ({ratios(straight)})",
                f"{curved[-1]:.1e} ({ratios(curved)})",
                ", ".join(f"{g['max_rel_area_error']:.1e}" for g in regions)
                + f" ({ratios([g['max_rel_area_error'] for g in regions])})",
                sum(g["non_simple_rings"] for g in regions),
                f"{fd['n_within_threshold']} within, {fd['n_crossing']} crossing, "
                f"Richardson max {fd['richardson_max_rel_error']:.1e}",
                f"{c['boundary']['max_offset_from_side']:.0e}",
            ]
        )

        bo, ao = b["optimization_cloud"], a["optimization_cloud"]
        quadrature.append(
            [
                name,
                f"{bo['nodes_per_cell_axis']}/{be['nodes_per_cell_axis']}",
                f"{bo['E0']:.6f} / {be['E0']:.6f}",
                f"{ao['E0']:.6f} / {ae['E0']:.6f}",
                f"{abs(ao['E0'] - ae['E0']) / ae['E0']:.1e}",
                f"{r['outcome']}, {r['n_accepted_steps']} steps",
                f"{r['seconds']['optimize']:.0f} s + {r['seconds']['replay']:.0f} s",
                r["settings"].get("threads", 1),
            ]
        )

    print("Effect sizes (evaluation cloud):\n")
    print(
        table(
            [
                "Run",
                "Modes",
                "E0 rectangle",
                "E0 land",
                "E0 ocean",
                "mean e rectangle",
                "mean e land",
                "mean e ocean",
                "e p99 land",
            ],
            effect,
        )
    )
    if weighting:
        print(
            "\nBoth weightings (evaluation cloud; the run minimized "
            "E_alpha):\n"
        )
        print(
            table(
                [
                    "Run",
                    "Modes",
                    "alpha",
                    "E0 rectangle",
                    "E1 rectangle",
                    "E0 land",
                    "E1 land",
                    "E0 ocean",
                    "E1 ocean",
                ],
                weighting,
            )
        )
    if subsets:
        print(
            "\nDistortion by subset (evaluation cloud; e = ln(a/b), "
            "ẽ = 2 arcsin((a − b)/(a + b))):\n"
        )
        print(
            table(
                [
                    "Run",
                    "Subset",
                    "Share",
                    "E0",
                    "mean e",
                    "mean ẽ",
                    "e p50",
                    "e p90",
                    "e p99",
                    "ẽ p99",
                ],
                subsets,
            )
        )
    print("\nCorrectness checks (refinement ratios in parentheses):\n")
    print(
        table(
            [
                "Run",
                "max det drift",
                "Mesh, finer cells",
                "Mesh, more samples/edge",
                "Region areas, δ = 0.5, 0.25, 0.125",
                "Non-simple rings",
                "FD Jacobian",
                "Boundary offset",
            ],
            checks,
        )
    )
    print("\nQuadrature and cost:\n")
    print(
        table(
            [
                "Run",
                "q opt/eval",
                "E0(T0) opt / eval",
                "E0(T) opt / eval",
                "rel. gap opt vs eval, after",
                "Outcome",
                "Optimize + replay",
                "Threads",
            ],
            quadrature,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
