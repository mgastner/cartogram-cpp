#!/usr/bin/env python3
"""Draw the explanatory figure of the paper as an Inkscape SVG.

Usage: paper_figure_concept.py PLOT_JSON CSV OUT_SVG [--steps STEPS_JSON]

PLOT_JSON is the *_conformality_plot.json of a run of the research binary on
the synthetic map in research/data/concept_figure (nine unit squares, with
the population densities in CSV). STEPS_JSON is the *_conformality_steps.json
of the same run. The figure has four panels:

  (a) the input map, with small circles (Tissot's indicatrix) on a lattice
      and one marked point x;
  (b) the cartogram T0, with the ellipses into which T0 turns the circles;
  (c) a velocity field of A, drawn on the arrow from (b) to (d). With
      --steps: streamlines of the first field of the run, over the regions
      of T0. Without: streamlines of one sine mode, and a small square
      carried by the flow;
  (d) the corrected map T = A o T0, with the ellipses of T.

The ellipse of a circle of radius r at x has the semi-axes r a and r b, with
a and b the singular values of the Jacobian matrix of the field cell that
contains x. The Jacobian matrix comes from the corners of the cell
(corners_t0 and corners_t). The lengths are in millimeters, and the layers
are Inkscape layers, so the figure can be edited in Inkscape.

The run for the paper (18 Nov 2026) used the research binary at eb3ca2ec,
with the settings of the paper runs, 32 x 32 modes, and ocean weight 0.1:

  cartogram nine_squares.geojson nine_squares_diag_1_3_9.csv -P 100000000 \\
    --conf_optimize --conf_alpha 0.5 --conf_lbfgs 8 --conf_integrator verlet \\
    --conf_substep_factor 1 --conf_eval_quad 4 --conf_modes 32 \\
    --conf_ocean_weight 0.1

The densities 1, 3, and 9 rise from the lower left to the upper right, so
ln g0 is not far from harmonic and A can remove part of the distortion. With
a dense square in the middle (a strict maximum of ln g0), E_1/2 on land fell
by only 1.5 % (ocean weight 1) or 3 % (0.1). With the gradient and ocean
weight 1, E_1/2 on land falls by 20 %, but the lattice ellipses hardly
change (mean ln(a/b) at the 64 lattice points +1 %); with ocean weight 0.1,
E_1/2 on land falls by 60 % (0.301 to 0.119) and the mean ln(a/b) at the
lattice points by 35 %, while E_1/2 on the ocean rises from 0.123 to 0.413.
"""

import argparse
import csv
import json
import math

import numpy as np
from contourpy import contour_generator

# Layout in millimeters (text width of the ICA template: 170 mm)
WIDTH = 170.0
HEIGHT = 53.0
PANEL = 41.0  # side of panels (a), (b), and (c)
SMALL = 25.0  # side of panel (d)
X_A, X_B, X_C = 1.0, 54.0, 127.0
X_D = X_B + PANEL + (X_C - X_B - PANEL - SMALL) / 2
Y_TOP = 6.0  # top of panels (a), (b), and (c)
Y_TOP_D = 6.0
ARROW_Y_T0 = Y_TOP + PANEL / 2
ARROW_Y_A = 40.0

FONT = "Times New Roman"
PT = 0.3528  # millimeters per point
TITLE_SIZE = 8 * PT
LABEL_SIZE = 7 * PT

# Colors: one blue ramp for the population density (light = sparse), a dark
# red for the ellipses, and neutral grays for lines
DENSITY_FILL = {1: "#eef3fa", 3: "#cddcef", 9: "#a9c3e4"}
ELLIPSE = "#b03a2e"
GRID = "#9a9a9a"
OUTLINE = "#333333"
FRAME = "#666666"

R_LATTICE = 6.0  # radius of the lattice circles, in lattice units
R_MARKED = 14.0  # radius of the circle at the marked point
MARKED_POINT = (165.0, 91.0)  # input position of the marked point
LATTICE_PER_SIDE = 8  # circles per side of the block of squares


def fmt(v):
    return f"{v:.3f}".rstrip("0").rstrip(".")


class Panel:
    """Maps lattice coordinates (y up) to millimeters (y down)."""

    def __init__(self, x0, y0, size, lx, ly):
        self.x0, self.y0 = x0, y0
        self.s = size / max(lx, ly)
        self.lx, self.ly = lx, ly

    def __call__(self, x, y):
        return self.x0 + x * self.s, self.y0 + (self.ly - y) * self.s


def path_d(points, close=False):
    parts = [
        ("M" if k == 0 else "L") + f"{fmt(x)},{fmt(y)}"
        for k, (x, y) in enumerate(points)
    ]
    return " ".join(parts) + (" Z" if close else "")


def text(x, y, spans, size, anchor="start", fill="#000000"):
    """spans: list of (string, style); style is a space-separated subset of
    i (italic), b (bold), sub (subscript), and op (operator font)."""
    out = [
        f'<text x="{fmt(x)}" y="{fmt(y)}" font-family="{FONT}" '
        f'font-size="{fmt(size)}" text-anchor="{anchor}" fill="{fill}" '
        'xml:space="preserve">'
    ]
    for s, style in spans:
        tokens = style.split()
        attrs = []
        if "op" in tokens:
            # The ring operator is missing from Times New Roman
            attrs.append('font-family="STIXGeneral"')
        if "b" in tokens:
            # Points and vectors are bold and upright, as in the paper
            attrs.append('font-weight="bold"')
        if "i" in tokens:
            attrs.append('font-style="italic"')
        if "sub" in tokens:
            attrs.append('baseline-shift="sub" font-size="70%"')
        out.append(f"<tspan {' '.join(attrs)}>{s}</tspan>")
    out.append("</text>")
    return "".join(out)


def arrow(x0, y, x1, stroke=0.3, head=1.6):
    return (
        f'<path d="M{fmt(x0)},{fmt(y)} L{fmt(x1 - head)},{fmt(y)}" '
        f'stroke="#000000" stroke-width="{stroke}" fill="none"/>'
        f'<path d="M{fmt(x1)},{fmt(y)} L{fmt(x1 - head)},{fmt(y - head / 2.6)} '
        f'L{fmt(x1 - head)},{fmt(y + head / 2.6)} Z" fill="#000000"/>'
    )


def layer(label, body):
    return (
        f'<g inkscape:groupmode="layer" inkscape:label="{label}" '
        f'id="layer-{label.replace(" ", "-")}">\n{body}\n</g>'
    )


def cell_jacobians(corners, nx, ny, cell):
    c = np.asarray(corners).reshape(ny + 1, nx + 1, 2)
    dx = ((c[:-1, 1:] - c[:-1, :-1]) + (c[1:, 1:] - c[1:, :-1])) / (2 * cell)
    dy = ((c[1:, :-1] - c[:-1, :-1]) + (c[1:, 1:] - c[:-1, 1:])) / (2 * cell)
    return c, np.stack([dx, dy], axis=-1)  # J[j, i] = [[dx0, dy0], [dx1, dy1]]


def image_of(c, x, y, cell):
    """Bilinear interpolation of the corner images at input position (x, y)."""
    u, w = x / cell, y / cell
    i, j = int(math.floor(u)), int(math.floor(w))
    fu, fw = u - i, w - j
    p = (
        (1 - fu) * (1 - fw) * c[j, i]
        + fu * (1 - fw) * c[j, i + 1]
        + (1 - fu) * fw * c[j + 1, i]
        + fu * fw * c[j + 1, i + 1]
    )
    return p, (j, i)


def ellipse_svg(panel, center, jac, r, stroke, width, fill="none", opacity=1.0):
    """Image of a circle of radius r under the linear map jac, at center."""
    t = np.linspace(0, 2 * np.pi, 73)
    circle = r * np.stack([np.cos(t), np.sin(t)])
    pts = (jac @ circle).T + center
    d = path_d([panel(x, y) for x, y in pts], close=True)
    return (
        f'<path d="{d}" stroke="{stroke}" stroke-width="{width}" '
        f'fill="{fill}" fill-opacity="{opacity}"/>'
    )


def semi_axes(jac):
    u, s, vt = np.linalg.svd(jac)
    return s[0], s[1], u  # columns of u: directions of the axes on the image


def first_field_panel(pd, steps, regions_t0, lx, ly):
    """Streamlines of the first velocity field of A, drawn over the regions of
    the cartogram T0, on which the field acts.

    The stream function is s = sum_mn c_mn sin(kx_m y1) sin(ky_n y2), with the
    coefficients of the first accepted step of the run (steps JSON, index
    (m - 1) * my + (n - 1)), kx_m = m pi / lx, ky_n = n pi / ly, and
    v = (ds/dy2, -ds/dy1). No square is carried here: the first step moves
    points by only a few lattice units (at most 10 of 256 in the paper run),
    so its image would cover the square."""
    mx, my = steps["modes"]
    if [float(lx), float(ly)] != [float(v) for v in steps["lattice"]]:
        raise SystemExit("ERROR: The steps JSON belongs to a different lattice.")
    first = steps["steps"][0]
    c = np.asarray(first["coeffs"], dtype=float).reshape(mx, my)
    kx = np.arange(1, mx + 1) * np.pi / lx
    ky = np.arange(1, my + 1) * np.pi / ly

    def vel(pts):
        pts = np.atleast_2d(pts)
        sx = np.sin(np.outer(pts[:, 0], kx))
        cx = np.cos(np.outer(pts[:, 0], kx)) * kx
        sy = np.sin(np.outer(pts[:, 1], ky))
        cy = np.cos(np.outer(pts[:, 1], ky)) * ky
        return np.stack([((sx @ c) * cy).sum(1), -((cx @ c) * sy).sum(1)], axis=1)

    body = [
        f'<path d="{path_d([pd(0, 0), pd(lx, 0), pd(lx, ly), pd(0, ly)], True)}" '
        f'fill="#ffffff" stroke="{FRAME}" stroke-width="0.2"/>'
    ]
    for poly in regions_t0:
        pts = [pd(x, y) for x, y in poly["ring"]]
        body.append(
            f'<path d="{path_d(pts, True)}" fill="#ececec" stroke="#b5b5b5" '
            'stroke-width="0.15" stroke-linejoin="round"/>'
        )

    g1 = np.linspace(0, lx, 257)
    g2 = np.linspace(0, ly, 257)
    S = np.sin(np.outer(g2, ky)) @ c.T @ np.sin(np.outer(kx, g1))
    U, W = np.meshgrid(g1, g2)
    gen = contour_generator(U, W, S)
    fractions = (0.15, 0.4, 0.65, 0.9)
    for extreme in (S.max(), S.min()):
        for frac in fractions:
            for line in gen.lines(frac * extreme):
                body.append(
                    f'<path d="{path_d([pd(x, y) for x, y in line[::2]], True)}" '
                    'stroke="#555555" stroke-width="0.2" fill="none"/>'
                )
        # Two arrowheads on the longest streamline of the second level
        line = max(gen.lines(fractions[1] * extreme), key=len)
        for k in (len(line) // 4, 3 * len(line) // 4):
            pnt = line[k]
            v = vel(pnt)[0]
            v = v / np.linalg.norm(v)
            size = 0.035 * max(lx, ly)
            tip, base = pnt + size * v, pnt - size * v
            nrm = np.array([-v[1], v[0]]) * 0.85 * size
            tri = [pd(*tip), pd(*(base + nrm)), pd(*(base - nrm))]
            body.append(f'<path d="{path_d(tri, True)}" fill="#555555"/>')

    return body


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("plot_json")
    ap.add_argument("csv")
    ap.add_argument("out_svg")
    ap.add_argument(
        "--steps",
        help="*_conformality_steps.json of the same run: draw the first "
        "velocity field of A in panel (c) instead of one sine mode",
    )
    args = ap.parse_args()

    d = json.load(open(args.plot_json))
    f = d["field"]
    nx, ny = f["shape"]
    cell = f["cell"]
    lx, ly = d["final_lattice"]
    with open(args.csv) as fh:
        density = {row["name"]: int(row["Population"]) for row in csv.DictReader(fh)}
    names = d["region_names"]

    c0, J0 = cell_jacobians(f["corners_t0"], nx, ny, cell)
    c1, J1 = cell_jacobians(f["corners_t"], nx, ny, cell)
    identity = np.eye(2)

    pa = Panel(X_A, Y_TOP, PANEL, lx, ly)
    pb = Panel(X_B, Y_TOP, PANEL, lx, ly)
    pc = Panel(X_C, Y_TOP, PANEL, lx, ly)

    # Block of squares on the input map
    orig = d["polygons"]["original"]
    xs = [p[0] for r in orig for p in r["ring"]]
    ys = [p[1] for r in orig for p in r["ring"]]
    bx0, bx1, by0, by1 = min(xs), max(xs), min(ys), max(ys)
    step = (bx1 - bx0) / LATTICE_PER_SIDE
    lattice = [
        (bx0 + (k + 0.5) * step, by0 + (l + 0.5) * step)
        for l in range(LATTICE_PER_SIDE)
        for k in range(LATTICE_PER_SIDE)
    ]

    layers = []

    # Regions, filled by population density
    body = []
    for key, panel in (("original", pa), ("t0", pb), ("t", pc)):
        for poly in d["polygons"][key]:
            fill = DENSITY_FILL[density[names[poly["region"]]]]
            pts = [panel(x, y) for x, y in poly["ring"]]
            body.append(f'<path d="{path_d(pts, True)}" fill="{fill}" stroke="none"/>')
    layers.append(layer("regions", "\n".join(body)))

    # Grid lines: every second line of the plot JSON (spacing 16)
    body = []
    spacing = d["grid_spacing"]
    for k in range(0, int(round(lx / spacing)) + 1, 2):
        v = k * spacing
        for a, b in (((v, 0), (v, ly)), ((0, v), (lx, v))):
            body.append(
                f'<path d="{path_d([pa(*a), pa(*b)])}" stroke="{GRID}" '
                'stroke-width="0.1" fill="none"/>'
            )
    for key, panel in (("t0", pb), ("t", pc)):
        lines = d["grid"][key]
        n_half = len(lines) // 2
        for idx, line in enumerate(lines):
            if (idx % n_half) % 2:
                continue
            pts = [panel(x, y) for x, y in line[::8]] + [panel(*line[-1])]
            body.append(
                f'<path d="{path_d(pts)}" stroke="{GRID}" stroke-width="0.1" fill="none"/>'
            )
    layers.append(layer("grid lines", "\n".join(body)))

    # Region outlines and frames
    body = []
    for key, panel in (("original", pa), ("t0", pb), ("t", pc)):
        for poly in d["polygons"][key]:
            pts = [panel(x, y) for x, y in poly["ring"]]
            body.append(
                f'<path d="{path_d(pts, True)}" fill="none" stroke="{OUTLINE}" '
                'stroke-width="0.3" stroke-linejoin="round"/>'
            )
        body.append(
            f'<path d="{path_d([panel(0, 0), panel(lx, 0), panel(lx, ly), panel(0, ly)], True)}" '
            f'fill="none" stroke="{FRAME}" stroke-width="0.2"/>'
        )
    layers.append(layer("outlines", "\n".join(body)))

    # Lattice of indicatrices
    body = []
    for x, y in lattice:
        body.append(ellipse_svg(pa, np.array([x, y]), identity * R_LATTICE, 1.0, ELLIPSE, 0.15))
        for c, J, panel in ((c0, J0, pb), (c1, J1, pc)):
            p, (j, i) = image_of(c, x, y, cell)
            body.append(ellipse_svg(panel, p, J[j, i], R_LATTICE, ELLIPSE, 0.15))
    layers.append(layer("indicatrices", "\n".join(body)))

    # Marked point x and its indicatrix
    mx, my = MARKED_POINT
    body = [
        ellipse_svg(pa, np.array([mx, my]), identity, R_MARKED, "#000000", 0.3, ELLIPSE, 0.35)
    ]
    px, py = pa(mx, my)
    body.append(f'<circle cx="{fmt(px)}" cy="{fmt(py)}" r="0.35" fill="#000000"/>')
    body.append(text(px - 2.4, py + 3.0, [("x", "b")], LABEL_SIZE, "end"))
    values = {}
    for tag, (c, J, panel) in (("t0", (c0, J0, pb)), ("t", (c1, J1, pc))):
        p, (j, i) = image_of(c, mx, my, cell)
        jac = J[j, i]
        a, b, u = semi_axes(jac)
        values[tag] = (a, b, a / b + b / a - 2, np.linalg.det(jac))
        body.append(ellipse_svg(panel, p, jac, R_MARKED, "#000000", 0.3, ELLIPSE, 0.35))
        for s_len, direction, label in ((a, u[:, 0], "a"), (b, u[:, 1], "b")):
            end = p + R_MARKED * s_len * direction
            q0, q1 = panel(*p), panel(*end)
            body.append(
                f'<path d="{path_d([q0, q1])}" stroke="#000000" stroke-width="0.2"/>'
            )
            if tag == "t0":
                off = np.array(q1) - np.array(q0)
                off = off / np.linalg.norm(off) * 1.4
                body.append(
                    text(q1[0] + off[0], q1[1] + off[1] + 0.9, [(label, "i")], LABEL_SIZE, "middle")
                )
        qx, qy = panel(*p)
        body.append(f'<circle cx="{fmt(qx)}" cy="{fmt(qy)}" r="0.3" fill="#000000"/>')
    layers.append(layer("marked point", "\n".join(body)))

    if args.steps:
        # Panel (d): the first velocity field of A, from the steps JSON
        pd = Panel(X_D, Y_TOP_D, SMALL, lx, ly)
        with open(args.steps) as fh:
            steps = json.load(fh)
        body = first_field_panel(pd, steps, d["polygons"]["t0"], lx, ly)
        layers.append(layer("flow panel", "\n".join(body)))
    else:
        # Panel (d): streamlines of one sine mode, s = sin(pi u) sin(pi w)
        q = 1.0
        pd = Panel(X_D, Y_TOP_D, SMALL, q, q)
        body = [
            f'<path d="{path_d([pd(0, 0), pd(q, 0), pd(q, q), pd(0, q)], True)}" '
            f'fill="#ffffff" stroke="{FRAME}" stroke-width="0.2"/>'
        ]
        n = 201
        g = np.linspace(0, 1, n)
        U, Wg = np.meshgrid(g, g)
        S = np.sin(np.pi * U) * np.sin(np.pi * Wg)
        gen = contour_generator(U, Wg, S)
        for level in (0.12, 0.35, 0.6, 0.85):
            for line in gen.lines(level):
                body.append(
                    f'<path d="{path_d([pd(x, y) for x, y in line[::2]], True)}" '
                    'stroke="#555555" stroke-width="0.2" fill="none"/>'
                )

        def vel(pnt):
            x, y = pnt
            return np.array(
                [np.pi * np.sin(np.pi * x) * np.cos(np.pi * y),
                 -np.pi * np.cos(np.pi * x) * np.sin(np.pi * y)]
            )

        # Arrowheads along the streamline of level 0.35
        for pnt in ((0.5, 0.115), (0.885, 0.5), (0.5, 0.885), (0.115, 0.5)):
            pnt = np.array(pnt)
            v = vel(pnt)
            v = v / np.linalg.norm(v)
            tip = pnt + 0.035 * v
            nrm = np.array([-v[1], v[0]])
            base = pnt - 0.035 * v
            tri = [pd(*tip), pd(*(base + 0.03 * nrm)), pd(*(base - 0.03 * nrm))]
            body.append(f'<path d="{path_d(tri, True)}" fill="#555555"/>')

        # A small square carried by the flow (RK4); the area stays the same
        def flow(pnt, tau, steps=400):
            h = tau / steps
            p = np.array(pnt, dtype=float)
            for _ in range(steps):
                k1 = vel(p)
                k2 = vel(p + h / 2 * k1)
                k3 = vel(p + h / 2 * k2)
                k4 = vel(p + h * k3)
                p = p + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            return p

        side = 0.13
        sx, sy = 0.24, 0.30
        edge = np.linspace(0, 1, 25)[:-1]
        square = (
            [(sx + side * t, sy) for t in edge]
            + [(sx + side, sy + side * t) for t in edge]
            + [(sx + side * (1 - t), sy + side) for t in edge]
            + [(sx, sy + side * (1 - t)) for t in edge]
        )
        tau = 0.12
        image = [flow(p, tau) for p in square]

        def area(pts):
            pts = np.asarray(pts)
            x, y = pts[:, 0], pts[:, 1]
            return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

        print(f"panel (d): square area {area(square):.6f}, image area {area(image):.6f}")
        for pts, fill in ((square, "#ffffff"), (image, ELLIPSE)):
            body.append(
                f'<path d="{path_d([pd(*p) for p in pts], True)}" stroke="#000000" '
                f'stroke-width="0.25" fill="{fill}" fill-opacity="0.35"/>'
            )
        layers.append(layer("flow panel", "\n".join(body)))

    # Arrows and labels
    a0, b0 = values["t0"][0], values["t0"][1]
    a1, b1 = values["t"][0], values["t"][1]
    m0, m1 = values["t0"][2], values["t"][2]
    g0 = values["t0"][3]
    print(
        f"marked point: g0 = {g0:.4f} (T: {values['t'][3]:.4f}); "
        f"T0: a = {a0:.3f}, b = {b0:.3f}, m = {m0:.3f}; "
        f"T: a = {a1:.3f}, b = {b1:.3f}, m = {m1:.3f}"
    )
    body = []
    body.append(arrow(X_A + PANEL + 1.5, ARROW_Y_T0, X_B - 1.5))
    body.append(text((X_A + PANEL + X_B) / 2, ARROW_Y_T0 - 1.3, [("T", "i"), ("0", "sub")], TITLE_SIZE, "middle"))
    body.append(arrow(X_B + PANEL + 1.5, ARROW_Y_A, X_C - 1.5))
    mid = (X_B + PANEL + X_C) / 2
    body.append(text(mid, ARROW_Y_A - 1.3, [("A", "i")], TITLE_SIZE, "middle"))
    body.append(
        text(mid, ARROW_Y_A + 3.6, [("det ", ""), ("DA", "i"), (" = 1", "")], LABEL_SIZE, "middle")
    )
    ty = Y_TOP - 2.0
    body.append(text(X_A, ty, [("(a) Input map", "")], TITLE_SIZE))
    body.append(
        text(X_B, ty, [("(b) Cartogram ", ""), ("T", "i"), ("0", "sub"), (", det ", ""),
                       ("DT", "i"), ("0", "sub"), (" = ", ""), ("g", "i"), ("0", "sub")],
             TITLE_SIZE)
    )
    body.append(
        text(X_C, ty, [("(d) ", ""), ("T", "i"), (" = ", ""), ("A", "i"), (" ∘ ", "op"),
                       ("T", "i"), ("0", "sub"), (", det ", ""), ("DT", "i"), (" = ", ""),
                       ("g", "i"), ("0", "sub")],
             TITLE_SIZE)
    )
    body.append(
        text(X_D + SMALL / 2, Y_TOP_D - 2.0, [("(c) First field " if args.steps else "(c) One field ", ""), ("v", "b"), (" of ", ""), ("A", "i")],
             TITLE_SIZE, "middle")
    )
    by = Y_TOP + PANEL + 3.6
    # Legend of the population density
    lx0 = X_A + 3.0
    body.append(text(lx0, by, [("Population density", "")], LABEL_SIZE))
    sx = lx0 + 24.0
    for dens in (1, 3, 9):
        body.append(
            f'<rect x="{fmt(sx)}" y="{fmt(by - 2.1)}" width="2.4" height="2.4" '
            f'fill="{DENSITY_FILL[dens]}" stroke="{OUTLINE}" stroke-width="0.15"/>'
        )
        body.append(text(sx + 3.0, by, [(str(dens), "")], LABEL_SIZE))
        sx += 5.4
    body.append(
        text(X_B + PANEL / 2, by,
             [("ab", "i"), (" = ", ""), ("g", "i"), ("0", "sub"), (f" = {g0:.2f},  ", ""),
              ("m", "i"), (f" = {m0:.2f}", "")],
             LABEL_SIZE, "middle")
    )
    body.append(
        text(X_C + PANEL / 2, by,
             [("ab", "i"), (" = ", ""), ("g", "i"), ("0", "sub"),
              (f" = {values['t'][3]:.2f},  ", ""), ("m", "i"), (f" = {m1:.2f}", "")],
             LABEL_SIZE, "middle")
    )
    layers.append(layer("labels", "\n".join(body)))

    svg = (
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
        'xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd" '
        f'width="{fmt(WIDTH)}mm" height="{fmt(HEIGHT)}mm" '
        f'viewBox="0 0 {fmt(WIDTH)} {fmt(HEIGHT)}" version="1.1">\n'
        '<sodipodi:namedview id="namedview" inkscape:document-units="mm" '
        'pagecolor="#ffffff" bordercolor="#666666"/>\n'
        "<title>Explanatory figure: T = A o T0</title>\n"
        + "\n".join(layers)
        + "\n</svg>\n"
    )
    with open(args.out_svg, "w") as fh:
        fh.write(svg)
    print(f"wrote {args.out_svg}")


if __name__ == "__main__":
    main()
