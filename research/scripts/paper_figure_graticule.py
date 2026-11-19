#!/usr/bin/env python3
"""Draw a four-panel map figure of the paper as an Inkscape SVG.

Usage: paper_figure_graticule.py PLOT_JSON OUT_SVG [--binned]

PLOT_JSON is the *_conformality_plot.json of a run of the research binary
on a map without insets. Panels, in a 2 x 2 grid at the text width:

  (a) the projection T0, with its graticule and the regions;
  (b) the corrected projection T = A o T0, with its graticule and the regions;
  (c) both graticules, T0 in gray and T in black;
  (d) the change m(T) - m(T0) of the excess MIPS energy, drawn on the input
      map, where E0 weights every point by its area.

The graticule is the image of the square grid of the input map (spacing in
the plot file, 8 lattice units). Every panel shows the regions plus a margin,
cut at the edge of the lattice. Panel (d) uses the ColorBrewer diverging
palette RdBu with 11 classes (blue: less distortion). Changes of m are
heavy-tailed, so the scale is symmetric-logarithmic: class breaks at
plus and minus five consecutive values of the sequence ..., 0.1, 0.3, 1, 3,
10, ..., the largest of them the largest such value below the 99th
percentile of |change|; the middle class around 0 is neutral. With --binned, every cell takes the color of its class; without,
the color varies continuously along the same scale.
"""

import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paper_figure_concept as concept  # noqa: E402  text, layer, path_d, fmt
import plot_distortion_field as pdf_  # noqa: E402  cells, OKLab ramp

WIDTH = 170.0
GAP_X = 6.0
GAP_Y = 9.0  # room for the titles of the second row
Y_TOP = 5.0
MARGIN = 0.06  # around the regions, relative to their extent
OUTLINE = "#333333"
FRAME = "#666666"
LAND = "#ededed"
GRID_ON_MAP = "#7a7a7a"
GRID_T0 = "#b3b3b3"
GRID_T = "#111111"
BAR_H = 1.8
# ColorBrewer RdBu, 11 classes, from blue (less distortion) to red
RDBU = ["#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7",
        "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"]


def breaks_for(abs_values, top=None):
    """Five positive class breaks for the change of m: consecutive values of
    the sequence ..., 0.1, 0.3, 1, 3, 10, ..., the largest of them the
    largest such value below the 99th percentile of |change|, or below top
    if given (to share one scale between figures)."""
    p99 = float(top) if top else float(np.percentile(abs_values, 99))
    k = math.floor(math.log10(p99))
    seq = [m * 10.0**j for j in range(k - 4, k + 1) for m in (1, 3)]
    seq = [v for v in seq if v <= p99 * (1 + 1e-12)]
    return seq[-5:]


def symlog_position(v, breaks):
    """Position in [0, 1] of v on the symmetric-logarithmic scale.

    The middle class [-c, c] maps linearly onto its 1/11 of the bar; each
    further class occupies another 1/11, with the value logarithmic inside.
    """
    c = breaks[0]
    a = abs(v)
    if a <= c:
        u = 0.5 * a / c
    else:
        edges = list(breaks) + [breaks[-1] * 3]
        u = 0.5
        for lo, hi in zip(edges[:-1], edges[1:]):
            if a >= hi:
                u += 1
                continue
            u += math.log(a / lo) / math.log(hi / lo)
            break
        u = min(u, 5.5)
    return 0.5 + math.copysign(u, v) / 11


class Box:
    def __init__(self, x0, y0, crop, width):
        (self.ax, self.ay), (bx, by) = crop
        self.s = width / (bx - self.ax)
        self.x0, self.y0, self.w = x0, y0, width
        self.h = (by - self.ay) * self.s

    def __call__(self, p):
        return (self.x0 + (p[0] - self.ax) * self.s, self.y0 + self.h - (p[1] - self.ay) * self.s)


def crop_of(rings, lattice):
    allv = np.vstack(rings)
    lo, hi = allv.min(axis=0), allv.max(axis=0)
    pad = MARGIN * max(hi - lo)
    return np.maximum(lo - pad, 0.0), np.minimum(hi + pad, np.asarray(lattice, dtype=float))


def clip_def(cid, box):
    pts = [(box.x0, box.y0), (box.x0 + box.w, box.y0), (box.x0 + box.w, box.y0 + box.h), (box.x0, box.y0 + box.h)]
    return f'<clipPath id="{cid}"><path d="{concept.path_d(pts, True)}"/></clipPath>'


def frame_svg(box):
    pts = [(box.x0, box.y0), (box.x0 + box.w, box.y0), (box.x0 + box.w, box.y0 + box.h), (box.x0, box.y0 + box.h)]
    return f'<path d="{concept.path_d(pts, True)}" fill="none" stroke="{FRAME}" stroke-width="0.2"/>'


def rings_svg(box, rings, fill, stroke, width):
    d = " ".join(concept.path_d([box(p) for p in r], True) for r in rings)
    return (f'<path d="{d}" fill="{fill}" fill-rule="evenodd" stroke="{stroke}" '
            f'stroke-width="{width}" stroke-linejoin="round"/>')


def grid_svg(box, lines, color, width):
    d = " ".join(concept.path_d([box(p) for p in line[::4]] + [box(line[-1])]) for line in lines)
    return f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}"/>'


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("plot_json")
    ap.add_argument("out_svg")
    ap.add_argument("--binned", action="store_true", help="color every cell by its class")
    ap.add_argument("--top", type=float,
                    help="largest class break of panel (d), to share one scale between figures")
    ap.add_argument("--grid-step", type=int, default=1,
                    help="draw every n-th line of the graticule (2 for dense maps such as the world)")
    args = ap.parse_args()

    d = json.load(open(args.plot_json))
    f = d["field"]
    initial, final = d["initial_lattice"], d["final_lattice"]
    polys = {k: [np.asarray(r["ring"], dtype=float) for r in d["polygons"][k]] for k in ("original", "t0", "t")}
    # The plot file lists the lines of one direction, then those of the other
    n_half = len(d["grid"]["t0"]) // 2
    grid = {k: [np.asarray(line, dtype=float) for idx, line in enumerate(d["grid"][k])
                if (idx % n_half) % args.grid_step == 0] for k in ("t0", "t")}

    crop_t = crop_of(polys["t0"] + polys["t"], final)
    crop_in = crop_of(polys["original"], initial)
    pw = (WIDTH - GAP_X) / 2
    top = [Box(0.0, Y_TOP, crop_t, pw), Box(pw + GAP_X, Y_TOP, crop_t, pw)]
    y2 = Y_TOP + top[0].h + GAP_Y
    bottom = [Box(0.0, y2, crop_t, pw), Box(pw + GAP_X, y2, crop_in, pw)]
    boxes = top + bottom

    # Change of m on the input map
    m0 = [None if v is None else float(v) for v in f["m_t0"]]
    m1 = [None if v is None else float(v) for v in f["m_t"]]
    delta = [None if a is None or b is None else b - a for a, b in zip(m0, m1)]
    quads = [np.asarray(q, dtype=float) for q in pdf_.input_quads(f, initial)]
    (cx0, cy0), (cx1, cy1) = crop_in
    shown = [(q, v) for q, v in zip(quads, delta)
             if v is not None and cx0 <= q[:, 0].mean() <= cx1 and cy0 <= q[:, 1].mean() <= cy1]
    breaks = breaks_for([abs(v) for _, v in shown], args.top)
    ramp = pdf_.Ramp(RDBU)

    def color(v):
        u = symlog_position(v, breaks)
        if args.binned:
            k = min(10, max(0, int(u * 11)))
            return RDBU[k]
        return ramp(u)

    defs = [clip_def(f"clip{k}", b) for k, b in enumerate(boxes)]
    layers = []

    body = []
    for k, (box, key) in enumerate(zip(boxes[:2], ("t0", "t"))):
        body.append(f'<g clip-path="url(#clip{k})">')
        body.append(rings_svg(box, polys[key], LAND, "none", 0))
        body.append(grid_svg(box, grid[key], GRID_ON_MAP, 0.1))
        body.append(rings_svg(box, polys[key], "none", OUTLINE, 0.18))
        body.append("</g>")
        body.append(frame_svg(box))
    layers.append(concept.layer("projections", "\n".join(body)))

    body = [f'<g clip-path="url(#clip2)">',
            grid_svg(boxes[2], grid["t0"], GRID_T0, 0.14),
            grid_svg(boxes[2], grid["t"], GRID_T, 0.12),
            "</g>", frame_svg(boxes[2])]
    layers.append(concept.layer("graticules", "\n".join(body)))

    bins = {}
    box = boxes[3]
    for q, v in shown:
        c = color(v)
        if not args.binned:
            # Quantize the continuous ramp finely, so that equal colors share a path
            c = ramp(round(symlog_position(v, breaks) * 110) / 110)
        bins.setdefault(c, []).append(concept.path_d([box(p) for p in q], True))
    body = ['<g clip-path="url(#clip3)">']
    for c, paths in bins.items():
        body.append(f'<path d="{" ".join(paths)}" fill="{c}" stroke="{c}" stroke-width="0.05" stroke-linejoin="round"/>')
    body.append(rings_svg(box, polys["original"], "none", OUTLINE, 0.15))
    body.append("</g>")
    body.append(frame_svg(box))
    layers.append(concept.layer("change of m", "\n".join(body)))

    # Titles
    body = []
    titles = [
        [("(a) Projection ", ""), ("T", "i"), ("0", "sub")],
        [("(b) ", ""), ("T", "i"), (" = ", ""), ("A", "i"), (" ∘ ", "op"), ("T", "i"), ("0", "sub")],
        [("(c) Graticules of ", ""), ("T", "i"), ("0", "sub"), (" (gray) and ", ""), ("T", "i"), (" (black)", "")],
        [("(d) Change of ", ""), ("m", "i"), (" on the input map", "")],
    ]
    for box, spans in zip(boxes, titles):
        body.append(concept.text(box.x0, box.y0 - 1.6, spans, concept.TITLE_SIZE))

    # Color bar under (d): 11 classes of equal width, labeled at the breaks
    bw = min(boxes[3].w, 72.0)
    bx = boxes[3].x0 + (boxes[3].w - bw) / 2
    by = max(boxes[2].y0 + boxes[2].h, boxes[3].y0 + boxes[3].h) + 2.0
    if args.binned:
        for k, c in enumerate(RDBU):
            body.append(f'<rect x="{concept.fmt(bx + k * bw / 11)}" y="{concept.fmt(by)}" '
                        f'width="{concept.fmt(bw / 11)}" height="{concept.fmt(BAR_H)}" fill="{c}"/>')
    else:
        stops = "".join(f'<stop offset="{concept.fmt(s)}" stop-color="{ramp(s)}"/>' for s in np.linspace(0, 1, 23))
        body.append(f'<defs><linearGradient id="dm-bar" x1="0" y1="0" x2="1" y2="0">{stops}</linearGradient></defs>'
                    f'<rect x="{concept.fmt(bx)}" y="{concept.fmt(by)}" width="{concept.fmt(bw)}" '
                    f'height="{concept.fmt(BAR_H)}" fill="url(#dm-bar)"/>')
    body.append(f'<rect x="{concept.fmt(bx)}" y="{concept.fmt(by)}" width="{concept.fmt(bw)}" '
                f'height="{concept.fmt(BAR_H)}" fill="none" stroke="{OUTLINE}" stroke-width="0.15"/>')
    ticks = [-b for b in reversed(breaks)] + [0.0] + list(breaks)
    for k, tv in enumerate(ticks):
        tx = bx + bw * symlog_position(tv, breaks) if tv else bx + bw / 2
        body.append(f'<path d="{concept.path_d([(tx, by + BAR_H), (tx, by + BAR_H + 0.6)])}" '
                    f'stroke="{OUTLINE}" stroke-width="0.15"/>')
        if k % 2 == 1 or tv == 0.0:
            # Label the powers of ten only, so that the labels fit; the
            # neutral middle class around 0 needs no label
            continue
        label = f"{tv:g}".replace("-", "−")
        body.append(concept.text(tx, by + BAR_H + 3.1, [(label, "")], concept.LABEL_SIZE * 0.9, "middle"))
    body.append(concept.text(
        bx + bw / 2, by + BAR_H + 6.4,
        [("m", "i"), ("(", ""), ("T", "i"), (") − ", ""), ("m", "i"), ("(", ""), ("T", "i"), ("0", "sub"),
         (") (blue: less distortion)", "")],
        concept.LABEL_SIZE, "middle"))
    layers.append(concept.layer("labels", "\n".join(body)))

    total_h = by + BAR_H + 8.0
    svg = (
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
        '<svg xmlns="http://www.w3.org/2000/svg" '
        'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
        'xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd" '
        f'width="{concept.fmt(WIDTH)}mm" height="{concept.fmt(total_h)}mm" '
        f'viewBox="0 0 {concept.fmt(WIDTH)} {concept.fmt(total_h)}" version="1.1">\n'
        '<sodipodi:namedview id="namedview" inkscape:document-units="mm" '
        'pagecolor="#ffffff" bordercolor="#666666"/>\n'
        "<title>Graticules of T0 and T</title>\n"
        f"<defs>{''.join(defs)}</defs>\n" + "\n".join(layers) + "\n</svg>\n"
    )
    with open(args.out_svg, "w") as fh:
        fh.write(svg)
    print(f"wrote {args.out_svg}: breaks of |change of m| {', '.join(f'{b:g}' for b in breaks)}"
          f"{' (binned)' if args.binned else ' (continuous)'}")


if __name__ == "__main__":
    main()
