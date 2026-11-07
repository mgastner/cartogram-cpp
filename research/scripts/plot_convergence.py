#!/usr/bin/env python3
"""Plot the convergence of conformality optimizer runs.

Reads the <inset>_conformality_opt.json file of each RUN_DIR and writes an
SVG with two panels that share the iteration axis:

  1. E0 relative to its value for T0 (optimization cloud),
  2. the gradient norm sqrt(R / R0) on a logarithmic axis.

Colour identifies the map (fixed order: USA, China, World, then others);
line style identifies the number of Fourier modes along the longer side.
If rsvg-convert is installed, a PNG and a PDF are written next to the SVG.
Only the Python standard library is used.

Usage: plot_convergence.py RUN_DIR [RUN_DIR ...] -o OUTPUT.svg
"""

import argparse
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

FONT = "Helvetica, Arial, sans-serif"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

# Categorical slots 1-3 of the reference palette (validated all-pairs)
MAP_COLOURS = [("USA", "#2a78d6"), ("China", "#eb6834"), ("World", "#1baf7a")]
DASHES = {16: "7,4", 32: "", 64: "1.5,3"}

PANEL_W = 360
PANEL_H = 220
MARGIN_L = 64
MARGIN_R = 20
MARGIN_T = 34
GAP = 70
MARGIN_B = 48
LEGEND_H = 30


def map_name(stem):
    for key, name in (("usa", "USA"), ("china", "China"), ("world", "World")):
        if key in stem:
            return name
    return stem.split("_")[0].capitalize()


def load(run_dir):
    files = sorted(Path(run_dir).glob("*_conformality_opt.json"))
    if len(files) != 1:
        sys.exit(f"{run_dir}: expected one *_conformality_opt.json")
    r = json.loads(files[0].read_text())
    history = r["history"]
    energy = [h["E"] for h in history] + [history[-1]["E_after"]]
    r0 = history[0]["R"]
    grad = [math.sqrt(h["R"] / r0) for h in history]
    return {
        "map": map_name(files[0].stem),
        "modes": r["modes"],
        "energy": [e / energy[0] for e in energy],
        "grad": grad,
    }


class Axes:
    def __init__(self, x0, y0, xmax, ylo, yhi, log):
        self.x0, self.y0 = x0, y0
        self.xmax, self.ylo, self.yhi, self.log = xmax, ylo, yhi, log

    def x(self, v):
        return self.x0 + v / self.xmax * PANEL_W

    def y(self, v):
        if self.log:
            t = (math.log10(v) - math.log10(self.ylo)) / (
                math.log10(self.yhi) - math.log10(self.ylo)
            )
        else:
            t = (v - self.ylo) / (self.yhi - self.ylo)
        return self.y0 + (1 - t) * PANEL_H


def text(x, y, s, size=10, anchor="middle", colour=INK_SECONDARY, extra=""):
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
        f'font-family="{FONT}" font-size="{size}" fill="{colour}" {extra}>{s}</text>'
    )


def x_ticks(xmax):
    step = 10 ** math.floor(math.log10(xmax / 4))
    for m in (1, 2, 2.5, 5, 10):
        if xmax / (m * step) <= 5:
            step *= m
            break
    return [i * step for i in range(int(xmax / step) + 1)]


def frame(ax, title, ylabel, yticks, fmt):
    out = [text(ax.x0, ax.y0 - 14, title, 12, "start", INK)]
    for t in yticks:
        y = ax.y(t)
        out.append(
            f'<line x1="{ax.x0:.1f}" y1="{y:.1f}" x2="{ax.x0 + PANEL_W:.1f}" '
            f'y2="{y:.1f}" stroke="{GRID}" stroke-width="0.6"/>'
        )
        out.append(text(ax.x0 - 6, y + 3.5, fmt(t), 10, "end"))
    for t in x_ticks(ax.xmax):
        x = ax.x(t)
        out.append(
            f'<line x1="{x:.1f}" y1="{ax.y0 + PANEL_H:.1f}" x2="{x:.1f}" '
            f'y2="{ax.y0 + PANEL_H + 3:.1f}" stroke="{AXIS}" stroke-width="0.8"/>'
        )
        out.append(text(x, ax.y0 + PANEL_H + 15, f"{t:g}"))
    out.append(
        f'<line x1="{ax.x0:.1f}" y1="{ax.y0 + PANEL_H:.1f}" '
        f'x2="{ax.x0 + PANEL_W:.1f}" y2="{ax.y0 + PANEL_H:.1f}" '
        f'stroke="{AXIS}" stroke-width="0.8"/>'
    )
    out.append(text(ax.x0 + PANEL_W / 2, ax.y0 + PANEL_H + 32, "Iteration", 11))
    cy = ax.y0 + PANEL_H / 2
    out.append(
        text(ax.x0 - 46, cy, ylabel, 11, "middle", INK_SECONDARY,
             f'transform="rotate(-90 {ax.x0 - 46:.1f} {cy:.1f})"')
    )
    return out


def polyline(ax, values, colour, dash):
    pts = " ".join(f"{ax.x(i):.1f},{ax.y(v):.1f}" for i, v in enumerate(values))
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<polyline points="{pts}" fill="none" stroke="{colour}" '
        f'stroke-width="2" stroke-linejoin="round" stroke-linecap="round"{dash_attr}/>'
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("run_dirs", nargs="+", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args()

    runs = [load(d) for d in args.run_dirs]
    order = {name: k for k, (name, _) in enumerate(MAP_COLOURS)}
    runs.sort(key=lambda r: (order.get(r["map"], 99), max(r["modes"])))
    colours = dict(MAP_COLOURS)
    extra = [m for m in dict.fromkeys(r["map"] for r in runs) if m not in colours]
    if extra:
        sys.exit(f"no colour slot for maps {extra}; at most three maps")

    xmax = max(len(r["energy"]) - 1 for r in runs)
    emin = min(min(r["energy"]) for r in runs)
    elo = math.floor(emin * 10) / 10
    glo = 10 ** math.floor(math.log10(min(min(r["grad"]) for r in runs)))

    width = MARGIN_L + PANEL_W + GAP + PANEL_W + MARGIN_R
    height = MARGIN_T + PANEL_H + MARGIN_B + LEGEND_H
    a = Axes(MARGIN_L, MARGIN_T, xmax, elo, 1.0, False)
    b = Axes(MARGIN_L + PANEL_W + GAP, MARGIN_T, xmax, glo, 1.0, True)

    body = []
    eticks = [round(elo + k * 0.1, 10) for k in range(int(round((1 - elo) / 0.1)) + 1)]
    body += frame(a, "E0 relative to T0", "E0 / E0(T0)", eticks, lambda t: f"{t:.1f}")
    gticks = [10.0 ** k for k in range(int(math.log10(glo)), 1)]
    body += frame(
        b, "Gradient norm", "√(R / R0)", gticks,
        lambda t: "1" if t >= 1 else f"10<tspan dy='-4' font-size='8'>{int(math.log10(t))}</tspan>",
    )
    for r in runs:
        colour = colours[r["map"]]
        dash = DASHES.get(max(r["modes"]), "3,2,1,2")
        body.append(polyline(a, r["energy"], colour, dash))
        body.append(polyline(b, r["grad"], colour, dash))

    # Legend: one entry per run, with a line sample in its colour and style
    x = MARGIN_L
    y = MARGIN_T + PANEL_H + MARGIN_B + 12
    for r in runs:
        label = f'{r["map"]}, {r["modes"][0]} × {r["modes"][1]} modes'
        dash = DASHES.get(max(r["modes"]), "3,2,1,2")
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        body.append(
            f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{x + 26:.1f}" y2="{y:.1f}" '
            f'stroke="{colours[r["map"]]}" stroke-width="2"{dash_attr}/>'
        )
        body.append(text(x + 32, y + 3.5, label, 10, "start", INK))
        x += 32 + 7 * len(label) + 18
        if x > width - 140:
            x = MARGIN_L
            y += 16
            height += 16

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" '
        f'height="{height:.0f}" viewBox="0 0 {width:.0f} {height:.0f}">\n'
        '<rect width="100%" height="100%" fill="#ffffff"/>\n'
        + "\n".join(body)
        + "\n</svg>\n"
    )
    args.output.write_text(svg)
    print(f"wrote {args.output}")
    converter = shutil.which("rsvg-convert")
    if converter:
        for fmt, flags in (("png", ["-z", "2"]), ("pdf", ["-f", "pdf"])):
            target = args.output.with_suffix("." + fmt)
            subprocess.run(
                [converter, *flags, "-o", str(target), str(args.output)], check=True
            )
            print(f"wrote {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
