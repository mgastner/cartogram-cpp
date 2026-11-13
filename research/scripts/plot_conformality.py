#!/usr/bin/env python3
"""Draw the square graticule and polygons under T0 and T = A o T0.

Reads the <inset>_conformality_plot.json file that cartogram writes with
--conformality (or --conf_optimize) and writes an SVG with four panels:

  1. input map with the square graticule,
  2. the cartogram T0 (fast flow-based algorithm),
  3. the corrected cartogram T = A o T0 (only if the optimizer ran),
  4. the graticule of T0 (grey) overlaid with that of T (black).

If rsvg-convert is installed, a PNG is written next to the SVG.
Only the Python standard library is used.

Usage: plot_conformality.py PLOT_JSON [-o OUTPUT.svg] [--width PIXELS]
"""

import argparse
import colorsys
import json
import shutil
import subprocess
import sys
from pathlib import Path

PANEL_GAP = 24
TITLE_HEIGHT = 28


def region_colour(index):
    """Pale, well-separated fill colour for a region index."""
    hue = (index * 0.6180339887) % 1.0
    r, g, b = colorsys.hls_to_rgb(hue, 0.82, 0.55)
    return f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"


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
        # SVG y grows downwards; lattice y grows upwards
        x = self.x0 + p[0] * self.scale
        y = self.y0 + (self.ly - p[1]) * self.scale
        return f"{x:.2f},{y:.2f}"

    def path(self, points, closed):
        if not points:
            return ""
        d = "M" + " L".join(self.point(p) for p in points)
        return d + (" Z" if closed else "")


def polygon_elements(panel, rings):
    """One filled path per region (even-odd fill keeps holes open)."""
    by_region = {}
    for item in rings:
        by_region.setdefault(item["region"], []).append(item["ring"])
    out = []
    for region, region_rings in sorted(by_region.items()):
        d = " ".join(panel.path(r, True) for r in region_rings)
        out.append(
            f'<path d="{d}" fill="{region_colour(region)}" '
            'fill-rule="evenodd" stroke="#555" stroke-width="0.5"/>'
        )
    return out


def line_elements(panel, lines, colour, width):
    return [
        f'<path d="{panel.path(line, False)}" fill="none" '
        f'stroke="{colour}" stroke-width="{width}"/>'
        for line in lines
    ]


def straight_graticule(lattice, spacing, samples_per_unit=1):
    """Lines of constant x and y in the input lattice."""
    lx, ly = lattice
    lines = []
    c = 0.0
    while c <= lx + 1e-9:
        lines.append([(min(c, lx), 0.0), (min(c, lx), ly)])
        c += spacing
    c = 0.0
    while c <= ly + 1e-9:
        lines.append([(0.0, min(c, ly)), (lx, min(c, ly))])
        c += spacing
    return lines


def frame(panel):
    return (
        f'<rect x="{panel.x0:.2f}" y="{panel.y0:.2f}" '
        f'width="{panel.width:.2f}" height="{panel.height:.2f}" '
        'fill="none" stroke="#000" stroke-width="0.8"/>'
    )


def title(panel, text):
    return (
        f'<text x="{panel.x0 + panel.width / 2:.2f}" '
        f'y="{panel.y0 - 9:.2f}" text-anchor="middle" '
        'font-family="Helvetica, Arial, sans-serif" font-size="15">'
        f"{text}</text>"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("plot_json", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument(
        "--width", type=float, default=520, help="panel width in pixels"
    )
    parser.add_argument(
        "--columns", type=int, default=2,
        help="panels per row (4 draws all panels in one row)",
    )
    args = parser.parse_args()

    data = json.loads(args.plot_json.read_text())
    output = args.output or args.plot_json.with_suffix(".svg")
    optimized = data["grid"]["t"] is not None
    initial = data["initial_lattice"]
    final = data["final_lattice"]
    spacing = data["grid_spacing"]

    panels_spec = [
        ("Input map", initial, "original"),
        ("T0: fast flow-based cartogram", final, "t0"),
    ]
    if optimized:
        panels_spec.append(("T = A ∘ T0: corrected cartogram", final, "t"))
        panels_spec.append(("Graticule: T0 (grey) and T (black)", final, "both"))

    width = args.width
    heights = [width * lat[1] / lat[0] for _, lat, _ in panels_spec]
    columns = max(1, args.columns)
    rows = (len(panels_spec) + columns - 1) // columns
    row_heights = [
        max(heights[r * columns : (r + 1) * columns]) for r in range(rows)
    ]
    total_w = columns * width + (columns + 1) * PANEL_GAP
    total_h = sum(h + TITLE_HEIGHT + PANEL_GAP for h in row_heights) + PANEL_GAP

    body = []
    y = PANEL_GAP
    for r in range(rows):
        for c in range(columns):
            k = r * columns + c
            if k >= len(panels_spec):
                break
            name, lattice, key = panels_spec[k]
            panel = Panel(PANEL_GAP + c * (width + PANEL_GAP), y, width, lattice)
            body.append(title(panel, name))
            if key == "original":
                body += polygon_elements(panel, data["polygons"]["original"])
                body += line_elements(
                    panel, straight_graticule(lattice, spacing), "#333", 0.5
                )
            elif key in ("t0", "t"):
                body += polygon_elements(panel, data["polygons"][key])
                body += line_elements(panel, data["grid"][key], "#333", 0.5)
            else:
                body += line_elements(panel, data["grid"]["t0"], "#aaa", 0.9)
                body += line_elements(panel, data["grid"]["t"], "#000", 0.5)
            body.append(frame(panel))
        y += row_heights[r] + TITLE_HEIGHT + PANEL_GAP

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{total_w:.0f}" '
        f'height="{total_h:.0f}" viewBox="0 0 {total_w:.0f} {total_h:.0f}">\n'
        f'<rect width="100%" height="100%" fill="#fff"/>\n'
        + "\n".join(body)
        + "\n</svg>\n"
    )
    output.write_text(svg)
    print(f"wrote {output}")

    converter = shutil.which("rsvg-convert")
    if converter:
        png = output.with_suffix(".png")
        subprocess.run(
            [converter, "-z", "2", "-o", str(png), str(output)], check=True
        )
        print(f"wrote {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
