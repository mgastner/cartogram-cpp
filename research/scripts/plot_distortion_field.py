#!/usr/bin/env python3
"""Draw the distortion field e = ln(a/b) of T0 and of T = A o T0.

Reads the <inset>_conformality_plot.json file that cartogram writes with
--conf_optimize and writes an SVG with three panels:

  1. e of T0, each field cell drawn at its image under T0,
  2. e of T, each field cell drawn at its image under T,
  3. the change e(T) - e(T0), drawn on the input map, where the objective
     E0 weights every point by its original area.

Panels 1 and 2 share one sequential scale (one hue, light to dark), capped at
the 99th percentile of the cell values of T0. Panel 3 uses a diverging scale
(blue: less distortion, red: more) capped symmetrically at the 99th
percentile of |change|. Caps are rounded up to round numbers, and values
beyond a cap take the end colour.

If rsvg-convert is installed, a PNG and a PDF are written next to the SVG.
Only the Python standard library is used.

Usage: plot_distortion_field.py PLOT_JSON [-o OUTPUT.svg] [--width PIXELS]
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
OUTLINE = "#52514e"
PANEL_GAP = 28
TITLE_HEIGHT = 26
BAR_HEIGHT = 10
BAR_BLOCK = 48

# Sequential blue ramp, steps 100 to 700 of the reference palette
SEQUENTIAL = [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
    "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
    "#0d366b",
]
# Diverging blue <-> red with a neutral grey midpoint
DIVERGING_LOW = "#256abf"
DIVERGING_MID = "#f0efec"
DIVERGING_HIGH = "#c0392f"


def hex_to_rgb(h):
    return tuple(int(h[i : i + 2], 16) / 255 for i in (1, 3, 5))


def rgb_to_hex(rgb):
    return "#" + "".join(f"{round(min(1, max(0, c)) * 255):02x}" for c in rgb)


def to_linear(c):
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def to_srgb(c):
    c = max(0.0, c)
    return 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055


def to_oklab(h):
    r, g, b = (to_linear(c) for c in hex_to_rgb(h))
    l = 0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b
    m = 0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b
    s = 0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b
    l, m, s = (math.copysign(abs(v) ** (1 / 3), v) for v in (l, m, s))
    return (
        0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
        1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
        0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s,
    )


def from_oklab(lab):
    L, a, b = lab
    l = (L + 0.3963377774 * a + 0.2158037573 * b) ** 3
    m = (L - 0.1055613458 * a - 0.0638541728 * b) ** 3
    s = (L - 0.0894841775 * a - 1.2914855480 * b) ** 3
    rgb = (
        4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
        -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
        -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s,
    )
    return rgb_to_hex(tuple(to_srgb(c) for c in rgb))


class Ramp:
    """Piecewise-linear interpolation between colour stops in OKLab."""

    def __init__(self, stops):
        self.stops = [to_oklab(h) for h in stops]

    def __call__(self, t):
        t = min(1.0, max(0.0, t)) * (len(self.stops) - 1)
        i = min(int(t), len(self.stops) - 2)
        f = t - i
        a, b = self.stops[i], self.stops[i + 1]
        return from_oklab(tuple(x + f * (y - x) for x, y in zip(a, b)))


def quantile(values, p):
    v = sorted(values)
    return v[min(len(v) - 1, int(p * (len(v) - 1) + 0.5))]


class Panel:
    """Maps lattice coordinates [0, lx] x [0, ly] to an SVG rectangle."""

    def __init__(self, x0, y0, width, lattice):
        self.lx, self.ly = lattice
        self.x0 = x0
        self.y0 = y0 + TITLE_HEIGHT
        self.scale = width / self.lx
        self.width = width
        self.height = self.ly * self.scale

    def point(self, p):
        x = self.x0 + p[0] * self.scale
        y = self.y0 + (self.ly - p[1]) * self.scale
        return f"{x:.2f},{y:.2f}"

    def path(self, points, closed=True):
        d = "M" + " L".join(self.point(p) for p in points)
        return d + (" Z" if closed else "")


def cell_quads(field, corners_key):
    """Corner images of every field cell, in the order of the cell values."""
    nx, ny = field["shape"]
    corners = field[corners_key]

    def corner(i, k):
        return corners[k * (nx + 1) + i]

    return [
        [corner(i, k), corner(i + 1, k), corner(i + 1, k + 1), corner(i, k + 1)]
        for k in range(ny)
        for i in range(nx)
    ]


def input_quads(field, lattice):
    nx, ny = field["shape"]
    cell = field["cell"]
    lx, ly = lattice
    quads = []
    for k in range(ny):
        for i in range(nx):
            x0, x1 = i * cell, min(lx, (i + 1) * cell)
            y0, y1 = k * cell, min(ly, (k + 1) * cell)
            quads.append([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    return quads


def field_elements(panel, quads, values, colour_of):
    # A stroke in the fill colour hides anti-aliasing seams between cells
    out = []
    for quad, v in zip(quads, values):
        if v is None:
            continue
        c = colour_of(v)
        out.append(
            f'<path d="{panel.path(quad)}" fill="{c}" stroke="{c}" '
            'stroke-width="0.35"/>'
        )
    return out


def outline_elements(panel, rings):
    d = " ".join(panel.path(item["ring"]) for item in rings)
    return [
        f'<path d="{d}" fill="none" stroke="{OUTLINE}" stroke-width="0.4" '
        'stroke-linejoin="round"/>'
    ]


def frame(panel):
    return (
        f'<rect x="{panel.x0:.2f}" y="{panel.y0:.2f}" width="{panel.width:.2f}" '
        f'height="{panel.height:.2f}" fill="none" stroke="{INK_SECONDARY}" '
        'stroke-width="0.6"/>'
    )


def title(panel, text):
    return (
        f'<text x="{panel.x0:.2f}" y="{panel.y0 - 9:.2f}" font-family="{FONT}" '
        f'font-size="13" fill="{INK}">{text}</text>'
    )


def colour_bar(x0, y0, width, colour_of, lo, hi, label, ticks, capped):
    n = 120
    out = []
    for j in range(n):
        v = lo + (hi - lo) * (j + 0.5) / n
        out.append(
            f'<rect x="{x0 + j * width / n:.2f}" y="{y0:.2f}" '
            f'width="{width / n + 0.3:.2f}" height="{BAR_HEIGHT}" '
            f'fill="{colour_of(v)}"/>'
        )
    for t in ticks:
        x = x0 + (t - lo) / (hi - lo) * width
        text = f"{t:g}".replace("-", "−")
        # End labels are anchored inward so that adjacent bars never collide
        anchor = "middle"
        if math.isclose(t, hi):
            anchor = "end"
            if capped:
                text = "≥ " + text
        elif math.isclose(t, lo):
            anchor = "start"
            if capped and lo < 0:
                text = "≤ " + text
        out.append(
            f'<line x1="{x:.2f}" y1="{y0 + BAR_HEIGHT:.2f}" x2="{x:.2f}" '
            f'y2="{y0 + BAR_HEIGHT + 3:.2f}" stroke="{INK_SECONDARY}" '
            'stroke-width="0.6"/>'
            f'<text x="{x:.2f}" y="{y0 + BAR_HEIGHT + 14:.2f}" '
            f'text-anchor="{anchor}" font-family="{FONT}" font-size="10" '
            f'fill="{INK_SECONDARY}">{text}</text>'
        )
    out.append(
        f'<text x="{x0 + width / 2:.2f}" y="{y0 + BAR_HEIGHT + 28:.2f}" '
        f'text-anchor="middle" font-family="{FONT}" font-size="11" '
        f'fill="{INK}">{label}</text>'
    )
    return out


def nice_ceiling(x):
    """Smallest value of the form m * 10^k >= x, m in a list of round numbers."""
    k = math.floor(math.log10(x))
    for m in (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10):
        if m * 10**k >= x * (1 - 1e-12):
            return m * 10**k
    return 10 ** (k + 1)


def nice_ticks(lo, hi, count=5):
    span = hi - lo
    step = 10 ** math.floor(math.log10(span / count))
    for m in (1, 2, 2.5, 5, 10):
        if span / (m * step) <= count:
            step *= m
            break
    first = math.ceil(lo / step - 1e-9) * step
    ticks = []
    t = first
    while t <= hi + 1e-9:
        ticks.append(round(t, 10))
        t += step
    return ticks


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("plot_json", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument(
        "--width", type=float, default=420, help="panel width in pixels"
    )
    args = parser.parse_args()

    data = json.loads(args.plot_json.read_text())
    field = data.get("field")
    if not field or field.get("e_t") is None:
        sys.exit("plot data has no distortion field after optimization")
    output = args.output or args.plot_json.with_name(
        args.plot_json.stem.replace("_conformality_plot", "") + "_distortion_field.svg"
    )
    initial = data["initial_lattice"]
    final = data["final_lattice"]

    e0 = field["e_t0"]
    e1 = field["e_t"]
    # Caps: 99th percentiles, rounded up so that the bar ends at a tick
    finite0 = [v for v in e0 if v is not None]
    cap = nice_ceiling(quantile(finite0, 0.99))
    delta = [
        None if a is None or b is None else b - a for a, b in zip(e0, e1)
    ]
    dcap = nice_ceiling(quantile([abs(v) for v in delta if v is not None], 0.99))

    sequential = Ramp(SEQUENTIAL)
    low_arm = Ramp([DIVERGING_MID, DIVERGING_LOW])
    high_arm = Ramp([DIVERGING_MID, DIVERGING_HIGH])

    def seq_colour(v):
        return sequential(v / cap)

    def div_colour(v):
        return high_arm(v / dcap) if v >= 0 else low_arm(-v / dcap)

    width = args.width
    panels = [
        ("e of T0, on T0", final),
        ("e of T = A ∘ T0, on T", final),
        ("e(T) − e(T0), on the input map", initial),
    ]
    height = max(width * lat[1] / lat[0] for _, lat in panels)
    total_w = len(panels) * width + (len(panels) + 1) * PANEL_GAP
    total_h = PANEL_GAP + TITLE_HEIGHT + height + 12 + BAR_BLOCK + PANEL_GAP

    body = []
    for k, (name, lattice) in enumerate(panels):
        panel = Panel(PANEL_GAP + k * (width + PANEL_GAP), PANEL_GAP, width, lattice)
        body.append(title(panel, name))
        if k == 0:
            body += field_elements(
                panel, cell_quads(field, "corners_t0"), e0, seq_colour
            )
            body += outline_elements(panel, data["polygons"]["t0"])
        elif k == 1:
            body += field_elements(
                panel, cell_quads(field, "corners_t"), e1, seq_colour
            )
            body += outline_elements(panel, data["polygons"]["t"])
        else:
            body += field_elements(
                panel, input_quads(field, initial), delta, div_colour
            )
            body += outline_elements(panel, data["polygons"]["original"])
        body.append(frame(panel))

    bar_y = PANEL_GAP + TITLE_HEIGHT + height + 12
    bar_w = 2 * width + PANEL_GAP
    body += colour_bar(
        PANEL_GAP,
        bar_y,
        bar_w,
        seq_colour,
        0.0,
        cap,
        "e = ln(a/b), mean over each cell (shared scale)",
        nice_ticks(0.0, cap),
        True,
    )
    dticks = [t for t in nice_ticks(-dcap, dcap, 4)]
    body += colour_bar(
        PANEL_GAP + 2 * (width + PANEL_GAP),
        bar_y,
        width,
        div_colour,
        -dcap,
        dcap,
        "change in e (blue: less distortion)",
        dticks,
        True,
    )

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{total_w:.0f}" '
        f'height="{total_h:.0f}" viewBox="0 0 {total_w:.0f} {total_h:.0f}">\n'
        '<rect width="100%" height="100%" fill="#ffffff"/>\n'
        + "\n".join(body)
        + "\n</svg>\n"
    )
    output.write_text(svg)
    print(f"wrote {output} (cap {cap:.3g}, change cap {dcap:.3g})")

    converter = shutil.which("rsvg-convert")
    if converter:
        for fmt, flags in (("png", ["-z", "2"]), ("pdf", ["-f", "pdf"])):
            target = output.with_suffix("." + fmt)
            subprocess.run(
                [converter, *flags, "-o", str(target), str(output)], check=True
            )
            print(f"wrote {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
