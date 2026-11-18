#!/usr/bin/env python3
"""Draw a map figure of the paper as an Inkscape SVG: the angular distortion
e = ln(a/b) of the projection T0 and of T = A o T0, and its change.

Usage: paper_figure_maps.py RUN_DIR OUT_SVG [--equal-area GEOJSON]

RUN_DIR holds the outputs of one run of the research binary: a
*_conformality_plot.json for every inset and the *_cartogram.geojson, which
holds T0 with the insets in place (it matches the T0 polygons of the plot
files to 1e-11). For a map with insets, --equal-area names the input map with
the insets in place, as cartogram-cpp writes it with -E. The figure has three
panels in a row, in the style of the explanatory figure (Times New Roman,
millimeters, Inkscape layers):

  (a) e of T0, every field cell drawn at its image under T0;
  (b) e of T, every field cell drawn at its image under T, same scale;
  (c) the change e(T) - e(T0), drawn on the input map, where E0 weights
      every point by its area on the input map.

Every inset is cropped to its regions with a small margin and drawn in its
own frame; the cells outside the regions are drawn too, because the
objective includes them. Each inset is placed by the similarity that maps
its regions onto their positions in cartogram-cpp's own layout: T0 and T
share the lattice of the inset, so one similarity serves panels (a) and (b).
Panels (a) and (b) share one sequential scale, capped at the 99th
percentile of e on T0; panel (c) uses a diverging scale (blue: less
distortion), capped symmetrically at the 99th percentile of |change|.
"""

import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paper_figure_concept as concept  # noqa: E402  text, layer, path_d, fmt
import plot_distortion_field as pdf_  # noqa: E402  ramps, cells, ticks

WIDTH = 170.0  # text width of the ICA template, in millimeters
GAP = 6.0  # between panels
Y_TOP = 5.0  # top of the panels
MARGIN = 0.04  # margin around the regions of an inset, relative to its size
BAR_H = 1.8
N_BINS = 48
TITLE_SIZE = concept.TITLE_SIZE
LABEL_SIZE = concept.LABEL_SIZE
OUTLINE = "#333333"
FRAME = "#666666"


def ring_vertices(geom):
    polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
    return [np.asarray(ring, dtype=float) for poly in polys for ring in poly]


def similarity(src, dst):
    """Scale s and translation t with dst = s src + t, from the bounding boxes."""
    sx = np.ptp(dst[:, 0]) / np.ptp(src[:, 0])
    sy = np.ptp(dst[:, 1]) / np.ptp(src[:, 1])
    if abs(sx / sy - 1) > 1e-6:
        sys.exit(f"ERROR: The inset is not a similarity image (scales {sx:.9g}, {sy:.9g}).")
    t = dst.min(axis=0) - sx * src.min(axis=0)
    return sx, t


def load_layout(path):
    """Vertices of every region in a GeoJSON, keyed by every property value."""
    features = json.load(open(path))["features"]
    by_value = {}
    for f in features:
        verts = np.vstack(ring_vertices(f["geometry"]))
        for v in f["properties"].values():
            if isinstance(v, str):
                by_value.setdefault(v, []).append(verts)
    return by_value


def fit(layout, names, rings):
    src = np.vstack([np.asarray(r["ring"], dtype=float) for r in rings])
    try:
        dst = np.vstack([v for n in names for v in layout[n]])
    except KeyError as err:
        sys.exit(f"ERROR: Region {err} is missing from the layout.")
    return similarity(src, dst)


class Inset:
    def __init__(self, plot_json):
        d = json.load(open(plot_json))
        self.name = d["inset"]
        self.field = d["field"]
        self.initial = d["initial_lattice"]
        self.final = d["final_lattice"]
        self.names = d["region_names"]
        self.polys = d["polygons"]
        self.e0 = self.field["e_t0"]
        self.e1 = self.field["e_t"]


def panel_items(inset, key, transform):
    """Cells and outlines of one inset for one panel, in layout coordinates."""
    s, t = transform
    if key == "original":
        quads = pdf_.input_quads(inset.field, inset.initial)
    else:
        quads = pdf_.cell_quads(inset.field, "corners_t0" if key == "t0" else "corners_t")
    quads = [np.asarray(q, dtype=float) * s + t for q in quads]
    rings = [np.asarray(r["ring"], dtype=float) * s + t for r in inset.polys[key]]
    allv = np.vstack(rings)
    lo, hi = allv.min(axis=0), allv.max(axis=0)
    pad = MARGIN * max(hi - lo)
    # The margin stops at the edge of the lattice, outside of which there are
    # no cells
    lattice = np.asarray(inset.initial if key == "original" else inset.final, dtype=float)
    crop = (np.maximum(lo - pad, t), np.minimum(hi + pad, lattice * s + t))
    return quads, rings, crop


class Box:
    """Maps layout coordinates (y up) into a panel rectangle in millimeters."""

    def __init__(self, x0, y0, extent, height):
        (self.ax, self.ay), (bx, by) = extent
        self.s = height / (by - self.ay)
        self.x0, self.y0, self.h = x0, y0, height
        self.w = (bx - self.ax) * self.s

    def __call__(self, p):
        return (self.x0 + (p[0] - self.ax) * self.s, self.y0 + self.h - (p[1] - self.ay) * self.s)


def cells_svg(box, quads, values, crop, ramp, lo, hi, clip_id):
    (cx0, cy0), (cx1, cy1) = crop
    bins = {}
    for q, v in zip(quads, values):
        if v is None:
            continue
        c = q.mean(axis=0)
        if not (cx0 <= c[0] <= cx1 and cy0 <= c[1] <= cy1):
            continue
        k = int(np.clip((v - lo) / (hi - lo) * N_BINS, 0, N_BINS - 1))
        bins.setdefault(k, []).append(concept.path_d([box(p) for p in q], True))
    out = [f'<g clip-path="url(#{clip_id})">']
    for k in sorted(bins):
        col = ramp((k + 0.5) / N_BINS)
        # A thin stroke in the same color hides the seams between cells
        out.append(
            f'<path d="{" ".join(bins[k])}" fill="{col}" stroke="{col}" '
            'stroke-width="0.05" stroke-linejoin="round"/>'
        )
    out.append("</g>")
    return out


def rect_points(crop):
    (x0, y0), (x1, y1) = crop
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def color_bar(x0, y0, width, ramp, lo, hi, label_spans, capped_high, capped_low=False):
    gid = f"bar{abs(hash((x0, y0))) % 10**8}"
    stops = "".join(
        f'<stop offset="{concept.fmt(s)}" stop-color="{ramp(s)}"/>' for s in np.linspace(0, 1, 17)
    )
    out = [
        f'<defs><linearGradient id="{gid}" x1="0" y1="0" x2="1" y2="0">{stops}</linearGradient></defs>',
        f'<rect x="{concept.fmt(x0)}" y="{concept.fmt(y0)}" width="{concept.fmt(width)}" '
        f'height="{concept.fmt(BAR_H)}" fill="url(#{gid})" stroke="{OUTLINE}" stroke-width="0.15"/>',
    ]
    for tv in pdf_.nice_ticks(lo, hi):
        tx = x0 + width * (tv - lo) / (hi - lo)
        out.append(
            f'<path d="{concept.path_d([(tx, y0 + BAR_H), (tx, y0 + BAR_H + 0.6)])}" '
            f'stroke="{OUTLINE}" stroke-width="0.15"/>'
        )
        txt = f"{tv:g}".replace("-", "−")
        if capped_high and abs(tv - hi) < 1e-9:
            txt = "≥ " + txt
        if capped_low and abs(tv - lo) < 1e-9:
            txt = "≤ " + txt
        out.append(concept.text(tx, y0 + BAR_H + 3.1, [(txt, "")], LABEL_SIZE, "middle"))
    out.append(concept.text(x0 + width / 2, y0 + BAR_H + 6.4, label_spans, LABEL_SIZE, "middle"))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("out_svg")
    ap.add_argument("--equal-area", help="input map with the insets in place (cartogram -E)")
    args = ap.parse_args()

    insets = [Inset(f) for f in sorted(glob.glob(os.path.join(args.run_dir, "*_conformality_plot.json")))]
    if not insets:
        sys.exit("ERROR: No *_conformality_plot.json in the run directory.")
    if len(insets) > 1:
        if not args.equal_area:
            sys.exit("ERROR: A map with insets needs --equal-area.")
        cart = load_layout(glob.glob(os.path.join(args.run_dir, "*_cartogram.geojson"))[0])
        ea = load_layout(args.equal_area)
        place_t = {i.name: fit(cart, i.names, i.polys["t0"]) for i in insets}
        place_in = {i.name: fit(ea, i.names, i.polys["original"]) for i in insets}
    else:
        place_t = place_in = {insets[0].name: (1.0, np.zeros(2))}

    # Scales from the cells that the figure shows
    panels = []
    for key, place in (("t0", place_t), ("t", place_t), ("original", place_in)):
        items = {i.name: panel_items(i, key, place[i.name]) for i in insets}
        lo = np.min([c[0] for _, _, c in items.values()], axis=0)
        hi = np.max([c[1] for _, _, c in items.values()], axis=0)
        panels.append((key, items, (lo, hi)))

    def shown(values, key):
        out = []
        for i in insets:
            quads, _, crop = dict((k, it) for k, it, _ in panels)[key][i.name]
            (cx0, cy0), (cx1, cy1) = crop
            for q, v in zip(quads, values(i)):
                c = q.mean(axis=0)
                if v is not None and cx0 <= c[0] <= cx1 and cy0 <= c[1] <= cy1:
                    out.append(v)
        return out

    e0_shown = shown(lambda i: i.e0, "t0")
    cap = pdf_.nice_ceiling(pdf_.quantile(e0_shown, 0.99))
    delta = lambda i: [None if a is None or b is None else b - a for a, b in zip(i.e0, i.e1)]
    d_shown = shown(delta, "original")
    dcap = pdf_.nice_ceiling(pdf_.quantile([abs(v) for v in d_shown], 0.99))
    seq = pdf_.Ramp(pdf_.SEQUENTIAL)
    div = pdf_.Ramp([pdf_.DIVERGING_LOW, pdf_.DIVERGING_MID, pdf_.DIVERGING_HIGH])

    # Common height: the three panels and two gaps fill the text width
    aspects = [(hi[0] - lo[0]) / (hi[1] - lo[1]) for _, _, (lo, hi) in panels]
    height = (WIDTH - 2 * GAP) / sum(aspects)
    x = 0.0
    boxes = []
    for (key, items, extent), a in zip(panels, aspects):
        boxes.append(Box(x, Y_TOP, extent, height))
        x += a * height + GAP

    layers = []
    defs = []
    for p, ((key, items, extent), box) in enumerate(zip(panels, boxes)):
        body = []
        for i in insets:
            quads, rings, crop = items[i.name]
            cid = f"clip-{key}-{i.name}"
            defs.append(
                f'<clipPath id="{cid}"><path d="{concept.path_d([box(q) for q in rect_points(crop)], True)}"/></clipPath>'
            )
            if key == "original":
                body += cells_svg(box, quads, delta(i), crop, div, -dcap, dcap, cid)
            else:
                body += cells_svg(box, quads, i.e0 if key == "t0" else i.e1, crop, seq, 0.0, cap, cid)
        layers.append(concept.layer(f"cells ({'abc'[p]})", "\n".join(body)))

    body = []
    for (key, items, extent), box in zip(panels, boxes):
        for i in insets:
            quads, rings, crop = items[i.name]
            cid = f"clip-{key}-{i.name}"
            d = " ".join(concept.path_d([box(q) for q in r], True) for r in rings)
            body.append(
                f'<g clip-path="url(#{cid})"><path d="{d}" fill="none" stroke="{OUTLINE}" '
                'stroke-width="0.12" stroke-linejoin="round"/></g>'
            )
            body.append(
                f'<path d="{concept.path_d([box(q) for q in rect_points(crop)], True)}" '
                f'fill="none" stroke="{FRAME}" stroke-width="0.2"/>'
            )
    layers.append(concept.layer("outlines and frames", "\n".join(body)))

    # Titles and color bars
    body = []
    ty = Y_TOP - 1.6
    titles = [
        [("(a) Projection ", ""), ("T", "i"), ("0", "sub")],
        [("(b) ", ""), ("T", "i"), (" = ", ""), ("A", "i"), (" ∘ ", "op"), ("T", "i"), ("0", "sub")],
        [("(c) Change of ", ""), ("e", "i"), (" on the input map", "")],
    ]
    for box, spans in zip(boxes, titles):
        body.append(concept.text(box.x0, ty, spans, TITLE_SIZE))
    by = Y_TOP + height + 2.0
    bx0, bx1 = boxes[0].x0, boxes[1].x0 + boxes[1].w
    bar_w = min(bx1 - bx0, 70.0)
    body += color_bar(
        (bx0 + bx1 - bar_w) / 2, by, bar_w, seq, 0.0, cap,
        [("e", "i"), (" = ln(", ""), ("a", "i"), ("/", ""), ("b", "i"), (")", "")],
        max(e0_shown) > cap,
    )
    cw = min(boxes[2].w, 45.0)
    body += color_bar(
        boxes[2].x0 + (boxes[2].w - cw) / 2, by, cw, div, -dcap, dcap,
        [("e", "i"), ("(", ""), ("T", "i"), (") − ", ""), ("e", "i"), ("(", ""), ("T", "i"), ("0", "sub"),
         (") (blue: less distortion)", "")],
        max(d_shown) > dcap, min(d_shown) < -dcap,
    )
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
        "<title>Distortion of T0 and T</title>\n"
        f"<defs>{''.join(defs)}</defs>\n" + "\n".join(layers) + "\n</svg>\n"
    )
    with open(args.out_svg, "w") as fh:
        fh.write(svg)
    print(
        f"wrote {args.out_svg}: {len(insets)} inset(s), panel height {height:.1f} mm, "
        f"e cap {cap:g} (99th percentile on T0), change cap {dcap:g}"
    )


if __name__ == "__main__":
    main()
