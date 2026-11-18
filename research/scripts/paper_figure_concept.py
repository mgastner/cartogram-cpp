#!/usr/bin/env python3
"""Draw the explanatory figure of the paper as an Inkscape SVG.

Usage: paper_figure_concept.py PLOT_JSON CSV OUT_SVG [--steps STEPS_JSON]

PLOT_JSON is the *_conformality_plot.json of a run of the research binary on
the synthetic map in research/data/concept_figure (nine unit squares, with
their populations in CSV). STEPS_JSON is the *_conformality_steps.json of
the same run. The figure has four panels of equal size, so they share one
scale:

  (a) the input regions, with small circles (Tissot's indicatrix) on a
      lattice and one marked point x;
  (b) the projection T0, with the ellipses into which T0 turns the circles;
  (c) the projection T = A o T0, with the ellipses of T; the arrow A runs
      from (b) to (c);
  (d) behind a thin rule, a velocity field of A, which explains how A is
      built. With --steps: streamlines of the first field of the run, over
      the regions of T0. Without: streamlines of one sine mode, and a small
      square carried by the flow.

The land is shaded by the density rho recovered from the Jacobian of T0,
rho = rho_bar det DT0, on a logarithmic color bar with the densities of the
input squares as ticks. Each field cell keeps its color in (b) and (c),
because A keeps det DT0. Populations are extensive, so they get no color.
Leader lines tie the values below (b) and (c) to the highlighted ellipse,
whose semi-axes a and b are labeled in both panels. Labels and thin lines on
the map have a white casing.

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
PANEL = 35.5  # side of all four panels, so they share one scale
GAP_T0 = 9.0  # between (a) and (b), with the arrow T0
GAP_A = 13.0  # between (b) and (c), with the arrow A
GAP_FIELD = 4.0  # between (c) and (d), with a thin rule: (d) explains A
X_A = 1.0
X_B = X_A + PANEL + GAP_T0
X_T = X_B + PANEL + GAP_A  # panel (c), the map T
X_FIELD = X_T + PANEL + GAP_FIELD  # panel (d), a velocity field of A
HALO = 0.6  # width of the white casing of labels and thin lines
Y_TOP = 6.0  # top of the panels
ARROW_Y = Y_TOP + PANEL / 2
LABEL_GAP = 3.6  # from the bottom of the panels to the baseline of the labels
LABEL_HALF_WIDTH = 12.0  # leaders start within this distance of the center
HEIGHT = Y_TOP + PANEL + 7.5  # room for the color bar and its tick labels

FONT = "Times New Roman"
PT = 0.3528  # millimeters per point
TITLE_SIZE = 8 * PT
LABEL_SIZE = 7 * PT

# Colors: one blue ramp for the density rho on a logarithmic scale (light =
# sparse), a dark red for the ellipses, and neutral grays for lines
RHO_RANGE = (0.7, 12.0)  # densities at the ends of the color bar
RHO_TICKS = (1, 3, 9)  # the population densities of the input squares
RHO_COLORS = ("#f5f8fc", "#86aad9")  # sparse and dense end, mixed in OKLab
RHO_BINS = 64  # color steps of the shading
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


def text(x, y, spans, size, anchor="start", fill="#000000", halo=False):
    """spans: list of (string, style); style is a space-separated subset of
    i (italic), b (bold), sub (subscript), and op (operator font). With
    halo, a white outline under the glyphs sets the text off the map."""
    casing = (
        f'stroke="#ffffff" stroke-width="{fmt(HALO)}" stroke-linejoin="round" '
        'paint-order="stroke" '
        if halo
        else ""
    )
    out = [
        f'<text x="{fmt(x)}" y="{fmt(y)}" font-family="{FONT}" '
        f'font-size="{fmt(size)}" text-anchor="{anchor}" fill="{fill}" {casing}'
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


def cased_line(p0, p1, width=0.2):
    """A black line on a white casing, legible on the grid and the ellipses."""
    d = path_d([p0, p1])
    return (
        f'<path d="{d}" stroke="#ffffff" stroke-width="{fmt(width + HALO)}" '
        'stroke-linecap="round" fill="none"/>'
        f'<path d="{d}" stroke="#000000" stroke-width="{fmt(width)}" fill="none"/>'
    )


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


def ring_area(ring):
    r = np.asarray(ring, dtype=float)
    x, y = r[:, 0], r[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def inside_rings(points, rings):
    """Even-odd test of the points (n x 2) against the rings."""
    x, y = points[:, 0], points[:, 1]
    result = np.zeros(len(points), dtype=bool)
    for ring in rings:
        r = np.asarray(ring, dtype=float)
        x0, y0 = r[:, 0], r[:, 1]
        x1, y1 = np.roll(x0, 1), np.roll(y0, 1)
        for k in range(len(r)):
            if y0[k] == y1[k]:
                continue  # a horizontal edge never crosses the ray
            crosses = (y0[k] > y) != (y1[k] > y)
            x_cross = x0[k] + (y - y0[k]) * (x1[k] - x0[k]) / (y1[k] - y0[k])
            result ^= crosses & (x < x_cross)
    return result


# OKLab (Ottosson 2020), so that the ramp changes evenly in lightness
def _oklab(hex_color):
    srgb = np.array([int(hex_color[k:k + 2], 16) / 255 for k in (1, 3, 5)])
    rgb = np.where(srgb <= 0.04045, srgb / 12.92, ((srgb + 0.055) / 1.055) ** 2.4)
    lms = np.array([[0.4122214708, 0.5363325363, 0.0514459929],
                    [0.2119034982, 0.6806995451, 0.1073969566],
                    [0.0883024619, 0.2817188376, 0.6299787005]]) @ rgb
    return np.array([[0.2104542553, 0.7936177850, -0.0040720468],
                     [1.9779984951, -2.4285922050, 0.4505937099],
                     [0.0259040371, 0.7827717662, -0.8086757660]]) @ np.cbrt(lms)


def _hex_from_oklab(lab):
    lms = np.array([[1.0, 0.3963377774, 0.2158037573],
                    [1.0, -0.1055613458, -0.0638541728],
                    [1.0, -0.0894841775, -1.2914855480]]) @ lab
    rgb = np.array([[4.0767416621, -3.3077115913, 0.2309699292],
                    [-1.2684380046, 2.6097574011, -0.3413193965],
                    [-0.0041960863, -0.7034186147, 1.7076147010]]) @ lms ** 3
    rgb = np.clip(rgb, 0.0, 1.0)
    srgb = np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * rgb ** (1 / 2.4) - 0.055)
    r, g, b = (int(round(255 * v)) for v in srgb)
    return f"#{r:02x}{g:02x}{b:02x}"


def density_color(t):
    """Color of the density ramp at position t in [0, 1] (0 = sparse)."""
    lo, hi = (_oklab(c) for c in RHO_COLORS)
    return _hex_from_oklab(lo + (hi - lo) * min(max(float(t), 0.0), 1.0))


def density_position(rho):
    """Position of the density rho on the logarithmic color scale."""
    lo, hi = RHO_RANGE
    return (np.log(rho) - math.log(lo)) / (math.log(hi) - math.log(lo))


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
    pc = Panel(X_T, Y_TOP, PANEL, lx, ly)

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

    # Density rho recovered from the Jacobian of T0: rho = rho_bar g0, with
    # rho_bar the mean density of the regions, which are unit squares. Each
    # field cell on land keeps its color in (b) and (c), because A keeps g0.
    areas = [ring_area(poly["ring"]) for poly in orig]
    unit_area = float(np.mean(areas))
    rho_bar = sum(density[names[poly["region"]]] for poly in orig) * unit_area / sum(areas)
    rho = rho_bar * np.linalg.det(J0)
    xs_c = (np.arange(nx) + 0.5) * cell
    ys_c = (np.arange(ny) + 0.5) * cell
    centers = np.stack(np.meshgrid(xs_c, ys_c), axis=-1).reshape(-1, 2)
    land = inside_rings(centers, [poly["ring"] for poly in orig]).reshape(ny, nx)
    bins = np.clip((density_position(rho) * RHO_BINS).astype(int), 0, RHO_BINS - 1)
    corners_in = np.stack(
        np.meshgrid(np.arange(nx + 1) * cell, np.arange(ny + 1) * cell), axis=-1
    )
    print(
        f"density: rho_bar = {rho_bar:.3f}; on land {rho[land].min():.3f} "
        f"to {rho[land].max():.3f}"
    )
    body = []
    for corners, panel in ((corners_in, pa), (c0, pb), (c1, pc)):
        for k in np.unique(bins[land]):
            color = density_color((k + 0.5) / RHO_BINS)
            quads = []
            for j, i in zip(*np.nonzero(land & (bins == k))):
                quad = (corners[j, i], corners[j, i + 1], corners[j + 1, i + 1], corners[j + 1, i])
                quads.append(path_d([panel(*q) for q in quad], True))
            # A thin stroke in the same color hides the seams between cells
            body.append(
                f'<path d="{" ".join(quads)}" fill="{color}" stroke="{color}" '
                'stroke-width="0.05" stroke-linejoin="round"/>'
            )
    layers.append(layer("density", "\n".join(body)))

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
    body.append(text(px - 2.4, py + 3.0, [("x", "b")], LABEL_SIZE, "end", halo=True))
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
            body.append(cased_line(q0, q1))
            # Labeled in (b) and (c): from (b) to (c), a shrinks, b grows,
            # and ab stays the same
            off = np.array(q1) - np.array(q0)
            off = off / np.linalg.norm(off) * 1.4
            body.append(
                text(q1[0] + off[0], q1[1] + off[1] + 0.9, [(label, "i")], LABEL_SIZE,
                     "middle", halo=True)
            )
        qx, qy = panel(*p)
        body.append(f'<circle cx="{fmt(qx)}" cy="{fmt(qy)}" r="0.3" fill="#000000"/>')
        # Leader from the values below the panel to the highlighted ellipse,
        # so the values read as those of this ellipse, not of the whole map
        left = panel.x0 + PANEL / 2 - LABEL_HALF_WIDTH
        right = panel.x0 + PANEL / 2 + LABEL_HALF_WIDTH
        start = np.array([min(max(qx, left), right), Y_TOP + PANEL + LABEL_GAP - 2.4])
        center = np.array([qx, qy])
        direction = (start - center) / np.linalg.norm(start - center)
        m_mm = panel.s * np.diag([1.0, -1.0]) @ (R_MARKED * jac)
        end = center + direction / np.linalg.norm(np.linalg.solve(m_mm, direction))
        body.append(cased_line(start, end, 0.25))
    layers.append(layer("marked point", "\n".join(body)))

    if args.steps:
        # Panel (d): the first velocity field of A, from the steps JSON
        pd = Panel(X_FIELD, Y_TOP, PANEL, lx, ly)
        with open(args.steps) as fh:
            steps = json.load(fh)
        body = first_field_panel(pd, steps, d["polygons"]["t0"], lx, ly)
        layers.append(layer("flow panel", "\n".join(body)))
    else:
        # Panel (d): streamlines of one sine mode, s = sin(pi u) sin(pi w)
        q = 1.0
        pd = Panel(X_FIELD, Y_TOP, PANEL, q, q)
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
    body.append(arrow(X_A + PANEL + 1.0, ARROW_Y, X_B - 1.0))
    body.append(text((X_A + PANEL + X_B) / 2, ARROW_Y - 1.3, [("T", "i"), ("0", "sub")], TITLE_SIZE, "middle"))
    # A maps (b) to (c); panel (d) explains how A is built, so it stands
    # apart, behind a thin rule
    body.append(arrow(X_B + PANEL + 1.0, ARROW_Y, X_T - 1.0))
    mid = (X_B + PANEL + X_T) / 2
    body.append(text(mid, ARROW_Y - 1.3, [("A", "i")], TITLE_SIZE, "middle"))
    body.append(
        text(mid, ARROW_Y + 3.6, [("det ", ""), ("DA", "i"), (" = 1", "")], LABEL_SIZE, "middle")
    )
    rule_x = X_T + PANEL + GAP_FIELD / 2
    body.append(
        f'<path d="{path_d([(rule_x, 1.0), (rule_x, HEIGHT - 1.0)])}" '
        f'stroke="{FRAME}" stroke-width="0.15" fill="none"/>'
    )
    ty = Y_TOP - 2.0
    body.append(text(X_A, ty, [("(a) Input regions", "")], TITLE_SIZE))
    body.append(
        text(X_B, ty, [("(b) Projection ", ""), ("T", "i"), ("0", "sub"), (", det ", ""),
                       ("DT", "i"), ("0", "sub"), (" = ", ""), ("g", "i"), ("0", "sub")],
             TITLE_SIZE)
    )
    body.append(
        text(X_T, ty, [("(c) ", ""), ("T", "i"), (" = ", ""), ("A", "i"), (" ∘ ", "op"),
                       ("T", "i"), ("0", "sub"), (", det ", ""), ("DT", "i"), (" = ", ""),
                       ("g", "i"), ("0", "sub")],
             TITLE_SIZE)
    )
    # With --steps, the field of the first of the steps whose flows make up A
    body.append(
        text(X_FIELD, ty, [("(d) Velocity field ", ""), ("v", "b"), (" of ", ""), ("A", "i")],
             TITLE_SIZE)
    )
    by = Y_TOP + PANEL + LABEL_GAP
    # Color bar of the density rho (logarithmic), with the densities of the
    # input squares as ticks
    lx0 = X_A + 1.0
    body.append(text(lx0, by, [("Density ", ""), ("ρ", "i")], LABEL_SIZE))
    bar_x0, bar_w, bar_top, bar_h = lx0 + 12.5, 21.0, by - 2.0, 1.8
    stops = "".join(
        f'<stop offset="{fmt(s)}" stop-color="{density_color(s)}"/>'
        for s in np.linspace(0.0, 1.0, 17)
    )
    body.append(
        f'<defs><linearGradient id="rho-ramp" x1="0" y1="0" x2="1" y2="0">{stops}'
        "</linearGradient></defs>"
        f'<rect x="{fmt(bar_x0)}" y="{fmt(bar_top)}" width="{fmt(bar_w)}" '
        f'height="{fmt(bar_h)}" fill="url(#rho-ramp)" stroke="{OUTLINE}" '
        'stroke-width="0.15"/>'
    )
    for tick in RHO_TICKS:
        tx = bar_x0 + bar_w * float(density_position(tick))
        body.append(
            f'<path d="{path_d([(tx, bar_top + bar_h), (tx, bar_top + bar_h + 0.6)])}" '
            f'stroke="{OUTLINE}" stroke-width="0.15"/>'
        )
        body.append(text(tx, bar_top + bar_h + 2.9, [(str(tick), "")], LABEL_SIZE, "middle"))
    body.append(
        text(X_B + PANEL / 2, by,
             [("ab", "i"), (" = ", ""), ("g", "i"), ("0", "sub"), (f" = {g0:.2f},  ", ""),
              ("m", "i"), (f" = {m0:.2f}", "")],
             LABEL_SIZE, "middle")
    )
    body.append(
        text(X_T + PANEL / 2, by,
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
