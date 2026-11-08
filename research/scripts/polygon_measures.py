#!/usr/bin/env python3
"""Polygon-based shape measures of the input map, T0 and T = A o T0.

Reads the <inset>_conformality_plot.json file that cartogram writes with
--conf_optimize and compares every region of a cartogram (T0 or T) with the
same region on the input map.

By default (--metrics miaji), it uses the three deformation measures of
Miaji et al. (2026, Cartogr. Geogr. Inf. Sci. 53, 547-562, section 5.1.3),
which follow Sun (2020a). Each region (all its polygons and holes) of the
input map and of the cartogram is first scaled to unit area and translated
so that its centroid is at the origin. Then

  symmetric difference  Area(R1 minus R2) + Area(R2 minus R1) of the whole
                        normalized regions;
  Frechet distance      between the outer rings of the region's largest
                        polygon on the two maps;
  Hausdorff distance    between the same rings, one-sided: the largest
                        distance from a point on the input ring to the
                        nearest point on the cartogram ring.

Each is averaged over regions. The Hausdorff distance is one-sided because
that is what Miaji et al. computed: their section 5.1.3 defines it as "the
maximum of the shortest distances from any point on one boundary to the
nearest point on the other boundary", and Boost.Geometry's
discrete_hausdorff_distance, which they used, takes the maximum over the
first range only. With the input ring first, this reproduces their Table S4
values for cartogram-cpp on the USA and China maps; the symmetric distance
does not. The symmetric distance is written to the JSON output as well.

Miaji et al. take the largest polygon on each map separately, because one of
their methods drops small polygons; here the rings correspond one to one, so
the polygon with the largest area on the input map is paired with its own
images. Boost.Geometry measures Frechet and Hausdorff distances between
vertices only. Here both rings are first resampled to a maximum segment
length (--segment, in units of the normalized region, default 0.01), and
Hausdorff distances are measured from vertices to segments, so that the
values approximate the distances between the continuous curves and do not
depend on how densely T0 and T are sampled; --segment 0 uses the vertices as
given. The rings of the input map, T0 and T start at corresponding points,
so the Frechet distance pairs the curves from those starting points.

Miaji et al. do not optimize over rotations. With --rotation, each of the
three measures is also computed with rotations of the normalized cartogram
region about its centroid permitted, minimized separately for each measure.
On a grid of 72 equally spaced angles,

  --rotation local   descends from the angle 0 to the nearest local minimum,
                     the rotation by which the region has turned;
  --rotation global  takes the best angle on the whole circle.

Either is then refined by golden-section search within one grid spacing.
For the Frechet distance, the grid is evaluated on rings resampled at four
times the --segment length, and the refinement at --segment. Both searches
can mistake other deformations for rotations. The global minimum turns
roughly centrally symmetric regions upside down (on the USA T0, California
by -179 degrees, Illinois by 146 degrees), and the local search turns
squashed near-rectangles by almost 90 degrees (Colorado, Connecticut).

  --rotation map     rotates each cartogram region back by its own rotation
                     under the map, so the three measures share one angle.

The rotation is the least-squares (orthogonal Procrustes) rotation of the
region's field cells (the 2 x 2 cells of the plot data) from the input map
to the cartogram. It uses the map's own point correspondence, so it cannot
mistake a symmetric shape for a rotated one. The angles in the JSON output
(degrees, positive counterclockwise) are the rotations applied to the
cartogram region.

With --metrics 2018 (or all), it also uses the three measures that Gastner,
Seguy and More (2018, Table 1) took from the cartogram literature:

  alpha  aspect ratio: longer over shorter side of the minimum-area bounding
         rectangle of the region (over all rotations), averaged over regions;
  delta  Hamming distance: area in exactly one of the input region and the
         cartogram region, divided by the sum of their areas, after the
         cartogram region is scaled to the input region's area; minimized
         over translations and summed over regions (each term is in [0, 1]);
  theta  relative position error: angle between the line through the
         centroids of two regions on the input map and the same line on the
         cartogram, averaged over all pairs of regions and divided by pi,
         so that theta is in [0, 1].

These definitions follow section 4 of the 2018 SI. It averages over
"polygons" without saying whether a region with islands is one shape. By
default (--unit region), all polygons of a region form one shape, as in
Miaji et al. (2026). With --unit polygon, every polygon (an outer ring minus
its holes) is a separate shape, so each island is scored on its own. The
rings of the input map, T0 and T are listed in the same order, which pairs
each polygon with its images.

The translation that minimizes the Hamming distance is found by Nelder-Mead
from the alignment of the centroids, so it is a local minimum.

Plot data written before polygon rings were stored at full precision round
coordinates to 1e-4, which can make nearly touching edges cross. Such rings
are repaired with shapely.make_valid and counted in the output; data at full
precision should need no repair.

Requires shapely (pip install shapely), which brings numpy; everything else
is the standard library.

Usage: polygon_measures.py PLOT_JSON [PLOT_JSON ...]
           [--metrics miaji|2018|all] [--unit region|polygon]
           [--segment LENGTH] [--rotation local|global|map]
           [--json OUT.json]
"""

import argparse
import json
import math
import sys
from pathlib import Path

try:
    import numpy as np
    import shapely
    from shapely.affinity import rotate, scale, translate
    from shapely.geometry import LineString, Polygon
    from shapely.ops import unary_union
except ImportError:
    sys.exit("polygon_measures.py needs shapely: pip install shapely")


def signed_area(ring):
    s = 0.0
    x0, y0 = ring[0]
    for (x1, y1), (x2, y2) in zip(ring[1:], ring[2:]):
        s += (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    return 0.5 * s


def region_geometries(rings, n_regions):
    """Union of the outer rings minus the union of the holes, per region.

    Outer rings are counterclockwise and holes clockwise, as in the CGAL
    polygons that cartogram writes."""
    outers = [[] for _ in range(n_regions)]
    holes = [[] for _ in range(n_regions)]
    for item in rings:
        ring = item["ring"]
        if len(ring) < 3:
            continue
        target = outers if signed_area(ring) > 0 else holes
        target[item["region"]].append(Polygon(ring))
    geoms = []
    n_repaired = 0
    for o, h in zip(outers, holes):
        parts = []
        for p in o + h:
            if not p.is_valid:
                n_repaired += 1
                p = shapely.make_valid(p)
            parts.append(p)
        outer = unary_union(parts[: len(o)]) if o else None
        hole = unary_union(parts[len(o) :]) if h else None
        if outer is None:
            geoms.append(None)
        elif hole is None:
            geoms.append(outer)
        else:
            geoms.append(outer.difference(hole))
    return geoms, n_repaired


def polygon_geometries(rings):
    """One shape per outer ring, minus the holes of its region that lie in
    it, in the order of the rings. Returns (region, shape) pairs."""
    outers = []
    holes = []
    for item in rings:
        ring = item["ring"]
        if len(ring) < 3:
            continue
        p = Polygon(ring)
        if signed_area(ring) > 0:
            outers.append([item["region"], p, []])
        else:
            holes.append((item["region"], p))
    n_repaired = 0
    for entry in outers:
        if not entry[1].is_valid:
            n_repaired += 1
            entry[1] = shapely.make_valid(entry[1])
    for region, h in holes:
        if not h.is_valid:
            n_repaired += 1
            h = shapely.make_valid(h)
        pt = h.representative_point()
        # The innermost outer ring of the same region that contains the hole
        hosts = [e for e in outers if e[0] == region and e[1].contains(pt)]
        if hosts:
            min(hosts, key=lambda e: e[1].area)[2].append(h)
    shapes = [
        (region, outer.difference(unary_union(hs)) if hs else outer)
        for region, outer, hs in outers
    ]
    return shapes, n_repaired


def outer_rings(rings):
    """The outer rings, as (region, coordinates), in the order that
    polygon_geometries uses. The coordinates keep their starting point."""
    return [
        (item["region"], item["ring"])
        for item in rings
        if len(item["ring"]) >= 3 and signed_area(item["ring"]) > 0
    ]


def shapes_of(data, unit):
    """Names and shapes of the units (polygons or regions) on each map."""
    names = data["region_names"]
    geoms = {}
    repaired = 0
    if unit == "region":
        for key in ("original", "t0", "t"):
            geoms[key], r = region_geometries(data["polygons"][key], len(names))
            repaired += r
        return list(names), geoms, repaired
    regions = None
    for key in ("original", "t0", "t"):
        pairs, r = polygon_geometries(data["polygons"][key])
        repaired += r
        if regions is None:
            regions = [region for region, _ in pairs]
        elif [region for region, _ in pairs] != regions:
            sys.exit("rings of the input map, T0 and T do not correspond")
        geoms[key] = [g for _, g in pairs]
    count = {}
    unit_names = []
    for region in regions:
        count[region] = count.get(region, 0) + 1
        unit_names.append(f"{names[region]}#{count[region]}")
    return unit_names, geoms, repaired


def hull_points(geom):
    hull = geom.convex_hull
    if hull.geom_type != "Polygon":
        return None
    return list(hull.exterior.coords)[:-1]


def aspect_ratio(geom):
    """Longer over shorter side of the minimum-area bounding rectangle."""
    pts = hull_points(geom)
    if not pts or len(pts) < 3:
        return None
    best_area = math.inf
    best_ratio = None
    n = len(pts)
    for i in range(n):
        (x1, y1), (x2, y2) = pts[i], pts[(i + 1) % n]
        ex, ey = x2 - x1, y2 - y1
        length = math.hypot(ex, ey)
        if length == 0.0:
            continue
        ux, uy = ex / length, ey / length
        us = [px * ux + py * uy for px, py in pts]
        vs = [-px * uy + py * ux for px, py in pts]
        w = max(us) - min(us)
        h = max(vs) - min(vs)
        if w * h < best_area:
            best_area = w * h
            best_ratio = max(w, h) / min(w, h) if min(w, h) > 0 else math.inf
    return best_ratio


def unit_area(geom):
    """The geometry scaled about its centroid to unit area."""
    f = 1.0 / math.sqrt(geom.area)
    c = geom.centroid
    return scale(geom, xfact=f, yfact=f, origin=c)


def nelder_mead(f, x0, step, max_evals=300, tol=1e-9):
    simplex = [list(x0), [x0[0] + step, x0[1]], [x0[0], x0[1] + step]]
    values = [f(p) for p in simplex]
    evals = 3
    while evals < max_evals:
        order = sorted(range(3), key=lambda i: values[i])
        simplex = [simplex[i] for i in order]
        values = [values[i] for i in order]
        size = max(math.dist(simplex[0], p) for p in simplex[1:])
        if values[2] - values[0] <= tol and size <= 1e-6:
            break
        centroid = [(simplex[0][k] + simplex[1][k]) / 2 for k in range(2)]
        worst = simplex[2]
        reflected = [2 * centroid[k] - worst[k] for k in range(2)]
        fr = f(reflected)
        evals += 1
        if fr < values[0]:
            expanded = [3 * centroid[k] - 2 * worst[k] for k in range(2)]
            fe = f(expanded)
            evals += 1
            simplex[2], values[2] = (expanded, fe) if fe < fr else (reflected, fr)
        elif fr < values[1]:
            simplex[2], values[2] = reflected, fr
        else:
            contracted = [(centroid[k] + worst[k]) / 2 for k in range(2)]
            fc = f(contracted)
            evals += 1
            if fc < values[2]:
                simplex[2], values[2] = contracted, fc
            else:
                for i in (1, 2):
                    simplex[i] = [
                        (simplex[0][k] + simplex[i][k]) / 2 for k in range(2)
                    ]
                    values[i] = f(simplex[i])
                    evals += 1
    i = min(range(3), key=lambda j: values[j])
    return simplex[i], values[i]


def hamming(original, cartogram):
    """Minimum over translations of the area in exactly one shape divided by
    the sum of the areas, with both shapes scaled to unit area."""
    p = unit_area(original)
    q = unit_area(cartogram)
    pc, qc = p.centroid, q.centroid
    q = translate(q, pc.x - qc.x, pc.y - qc.y)

    def f(t):
        # (2 - 2 |p & q|) / (|p| + |q|) with |p| = |q| = 1
        return 1.0 - p.intersection(translate(q, t[0], t[1])).area

    best, value = nelder_mead(f, [0.0, 0.0], 0.05)
    # A restart from the best point guards against a collapsed simplex
    best, value = nelder_mead(f, best, 0.01)
    return value


def position_error(c_in, c_out):
    """Mean angle between centroid connections over all pairs, divided by
    pi."""
    total = 0.0
    n = 0
    for i in range(len(c_in)):
        for j in range(i + 1, len(c_in)):
            ax, ay = c_in[j][0] - c_in[i][0], c_in[j][1] - c_in[i][1]
            bx, by = c_out[j][0] - c_out[i][0], c_out[j][1] - c_out[i][1]
            total += abs(math.atan2(ax * by - ay * bx, ax * bx + ay * by))
            n += 1
    return total / (n * math.pi) if n else 0.0


def discrete_frechet(a, b):
    """Discrete Frechet distance between two vertex sequences (Eiter and
    Mannila 1994), computed along anti-diagonals of the coupling table.
    GEOS computes the same value but is too slow for rings with thousands of
    vertices."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n, m = len(a), len(b)
    prev2 = prev = None
    for k in range(n + m - 1):
        i = np.arange(max(0, k - m + 1), min(n - 1, k) + 1)
        j = k - i
        d = np.hypot(a[i, 0] - b[j, 0], a[i, 1] - b[j, 1])
        cur = np.full(n, np.inf)
        if k == 0:
            cur[0] = d[0]
        else:
            up = np.where(i >= 1, prev[np.maximum(i - 1, 0)], np.inf)
            left = np.where(j >= 1, prev[i], np.inf)
            diag = np.full(len(i), np.inf)
            if prev2 is not None:
                diag = np.where((i >= 1) & (j >= 1), prev2[np.maximum(i - 1, 0)], np.inf)
            cur[i] = np.maximum(d, np.minimum(np.minimum(up, left), diag))
        prev2, prev = prev, cur
    return float(prev[n - 1])


def one_sided_hausdorff(a, b):
    """Largest distance from a vertex of line a to line b."""
    return float(shapely.distance(shapely.points(np.asarray(a.coords)), b).max())


def golden_min(f, lo, hi, tol=1e-3):
    """Local minimum of f on [lo, hi] by golden-section search; returns the
    argument and the value."""
    g = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    while b - a > tol:
        if fc <= fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
    return (c, fc) if fc <= fd else (d, fd)


def min_over_rotation(f, coarse=None, n_grid=72, mode="local"):
    """Minimum of f(theta) over rotation angles theta (radians) and the angle
    in (-pi, pi]. On n_grid equally spaced angles (evaluated with coarse, a
    cheaper approximation of f, if given), mode "local" descends from the
    angle 0 to the nearest local minimum and mode "global" takes the best
    angle; either is refined by golden-section search within one grid
    spacing. The angle 0 is always a candidate, so the result never exceeds
    f(0)."""
    step = 2.0 * math.pi / n_grid
    g = coarse or f
    cache = {}

    def grid(k):
        k %= n_grid
        if k not in cache:
            cache[k] = g(k * step)
        return cache[k]

    if mode == "global":
        k = min(range(n_grid), key=grid)
    else:
        k = 0
        for _ in range(n_grid // 2):
            left, right = grid(k - 1), grid(k + 1)
            if min(left, right) >= grid(k):
                break
            k = k - 1 if left < right else k + 1
    theta, value = golden_min(f, (k - 1) * step, (k + 1) * step)
    f0 = f(0.0)
    if f0 <= value:
        return f0, 0.0
    return value, math.remainder(theta, 2.0 * math.pi)


def rotated_coords(coords, theta):
    """Coordinates rotated counterclockwise by theta about the origin."""
    xy = np.asarray(coords, dtype=float)
    c, s = math.cos(theta), math.sin(theta)
    return np.column_stack((c * xy[:, 0] - s * xy[:, 1], s * xy[:, 0] + c * xy[:, 1]))


def resampled(line, spacing):
    """Points at equal arc-length spacing (at most `spacing`) along line,
    including both ends."""
    n = max(2, int(math.ceil(line.length / spacing)) + 1)
    pts = shapely.line_interpolate_point(line, np.linspace(0.0, line.length, n))
    return shapely.get_coordinates(pts)


def rotation_measures(p, q, a, b, segment, mode="local", angle=None):
    """Symmetric difference, Frechet and one-sided Hausdorff distance after a
    rotation of the cartogram geometry (q, b) about the origin: each
    minimized over rotations (mode "local" or "global"), or all three at the
    given angle (mode "map", angle in radians). p, q: normalized regions;
    a, b: normalized outer rings (input, cartogram), already resampled at
    `segment`. Returns {measure: (value, angle in degrees)}, where the angle
    is the rotation applied to the cartogram geometry."""
    a_xy = np.asarray(a.coords)
    b_xy = np.asarray(b.coords)
    a_pts = shapely.points(a_xy)
    funcs = {
        "symdiff": lambda th: p.symmetric_difference(
            rotate(q, th, origin=(0, 0), use_radians=True)
        ).area,
        "frechet": lambda th: discrete_frechet(a_xy, rotated_coords(b_xy, th)),
        "hausdorff": lambda th: float(
            shapely.distance(a_pts, LineString(rotated_coords(b_xy, th))).max()
        ),
    }
    if mode == "map":
        return {m: (f(angle), math.degrees(angle)) for m, f in funcs.items()}
    coarse_spacing = 4.0 * segment if segment > 0 else 0.04
    a_coarse = resampled(a, coarse_spacing)
    b_coarse = resampled(b, coarse_spacing)
    coarse = {
        "frechet": lambda th: discrete_frechet(
            a_coarse, rotated_coords(b_coarse, th)
        )
    }
    out = {}
    for m, f in funcs.items():
        v, t = min_over_rotation(f, coarse=coarse.get(m), mode=mode)
        out[m] = (v, math.degrees(t))
    return out


def map_rotations(data, regions_original):
    """Rotation of each region under T0 and under T, from the map itself:
    the least-squares (orthogonal Procrustes) rotation that takes the
    centred centres of the field cells inside the region on the input map to
    their centred images, i.e. atan2(sum p x q, sum p . q). The images of
    the cell centres are the means of the images of the cell corners.
    Returns {key: [angle in radians or None]} for key in ("t0", "t"); None
    for a region with fewer than three cell centres."""
    f = data["field"]
    nx, ny = f["shape"]
    cell = f["cell"]
    jj, ii = np.mgrid[0:ny, 0:nx]
    xc = ((ii + 0.5) * cell).ravel()
    yc = ((jj + 0.5) * cell).ravel()
    masks = [
        None if g is None else shapely.contains_xy(g, xc, yc)
        for g in regions_original
    ]
    out = {}
    for key in ("t0", "t"):
        c = np.asarray(f[f"corners_{key}"], dtype=float).reshape(ny + 1, nx + 1, 2)
        img = 0.25 * (c[:-1, :-1] + c[:-1, 1:] + c[1:, :-1] + c[1:, 1:])
        img = img.reshape(-1, 2)
        angles = []
        for mask in masks:
            if mask is None or mask.sum() < 3:
                angles.append(None)
                continue
            px, py = xc[mask] - xc[mask].mean(), yc[mask] - yc[mask].mean()
            qx = img[mask, 0] - img[mask, 0].mean()
            qy = img[mask, 1] - img[mask, 1].mean()
            angles.append(
                math.atan2(np.sum(px * qy - py * qx), np.sum(px * qx + py * qy))
            )
        out[key] = angles
    return out


def normalized(geom, region):
    """geom under the map that takes region to unit area with its centroid
    at the origin."""
    f = 1.0 / math.sqrt(region.area)
    c = region.centroid
    return scale(translate(geom, -c.x, -c.y), xfact=f, yfact=f, origin=(0, 0))


def miaji_measures(data, segment, rotation=None):
    """Symmetric difference, Frechet and Hausdorff distance per region, after
    Miaji et al. (2026, section 5.1.3)."""
    names = data["region_names"]
    regions = {}
    rings = {}
    repaired = 0
    for key in ("original", "t0", "t"):
        regions[key], r = region_geometries(data["polygons"][key], len(names))
        repaired += r
        rings[key] = outer_rings(data["polygons"][key])
    order = [region for region, _ in rings["original"]]
    for key in ("t0", "t"):
        if [region for region, _ in rings[key]] != order:
            sys.exit("rings of the input map, T0 and T do not correspond")
    # The polygon of each region with the largest area on the input map
    largest = {}
    for k, (region, ring) in enumerate(rings["original"]):
        a = abs(signed_area(ring))
        if region not in largest or a > largest[region][1]:
            largest[region] = (k, a)
    turns = None
    if rotation == "map":
        if "field" not in data:
            sys.exit("--rotation map needs the distortion field in the plot data")
        turns = map_rotations(data, regions["original"])
    units = []
    for i, name in enumerate(names):
        if i not in largest or any(
            regions[key][i] is None or regions[key][i].area <= 0 for key in regions
        ):
            continue
        k = largest[i][0]
        ref = regions["original"][i]
        a = normalized(LineString(rings["original"][k][1]), ref)
        p = normalized(ref, ref)
        if segment > 0:
            a = shapely.segmentize(a, segment)
        entry = {"id": name}
        for key in ("t0", "t"):
            reg = regions[key][i]
            q = normalized(reg, reg)
            b = normalized(LineString(rings[key][k][1]), reg)
            if segment > 0:
                b = shapely.segmentize(b, segment)
            entry[f"symdiff_{key}"] = p.symmetric_difference(q).area
            entry[f"frechet_{key}"] = discrete_frechet(a.coords, b.coords)
            h_in = one_sided_hausdorff(a, b)
            entry[f"hausdorff_{key}"] = h_in
            entry[f"hausdorff_sym_{key}"] = max(h_in, one_sided_hausdorff(b, a))
            # Distance between the normalized starting points, to check that
            # the rings start at corresponding points
            entry[f"start_{key}"] = math.dist(a.coords[0], b.coords[0])
            if rotation:
                # In mode "map", undo the region's own rotation under the map
                # (no rotation if the region holds too few field cells)
                angle = None
                if turns is not None:
                    angle = -(turns[key][i] or 0.0)
                for m, (v, angle) in rotation_measures(
                    p, q, a, b, segment, rotation, angle
                ).items():
                    entry[f"{m}_rot_{key}"] = v
                    entry[f"{m}_angle_{key}"] = angle
        units.append(entry)
    names_avg = ["symdiff", "frechet", "hausdorff", "hausdorff_sym"]
    if rotation:
        names_avg += ["symdiff_rot", "frechet_rot", "hausdorff_rot"]
    maps = {}
    for key in ("t0", "t"):
        maps[key] = {
            m: sum(u[f"{m}_{key}"] for u in units) / len(units) for m in names_avg
        }
        if rotation:
            for m in ("symdiff", "frechet", "hausdorff"):
                angles = sorted(abs(u[f"{m}_angle_{key}"]) for u in units)
                maps[key][f"{m}_median_abs_angle"] = angles[len(angles) // 2]
    return {"n_regions": len(units), "maps": maps, "units": units}, repaired


def measures(path, unit, metrics="2018", segment=0.01, rotation=None):
    data = json.loads(Path(path).read_text())
    if data["polygons"]["t"] is None:
        sys.exit(f"{path}: no corrected polygons (run with --conf_optimize)")
    result = {"file": str(path)}
    if metrics in ("miaji", "all"):
        result["miaji"], repaired = miaji_measures(data, segment, rotation)
        result["miaji"]["segment"] = segment
        result["n_rings_repaired"] = repaired
        if metrics == "miaji":
            return result
    names, geoms, repaired = shapes_of(data, unit)
    n = len(names)
    keep = [
        i
        for i in range(n)
        if all(geoms[k][i] is not None and geoms[k][i].area > 0 for k in geoms)
    ]
    result.update(
        {
            "unit": unit,
            "n_units": len(keep),
            "n_rings_repaired": repaired,
            "maps": {},
            "units": [],
        }
    )
    centroids = {
        k: [(geoms[k][i].centroid.x, geoms[k][i].centroid.y) for i in keep]
        for k in geoms
    }
    per_unit = {i: {"id": names[i]} for i in keep}
    for key in ("original", "t0", "t"):
        ratios = []
        for i in keep:
            a = aspect_ratio(geoms[key][i])
            per_unit[i][f"aspect_{key}"] = a
            ratios.append(a)
        entry = {"alpha_mean": sum(ratios) / len(ratios)}
        if key != "original":
            deltas = []
            for i in keep:
                d = hamming(geoms["original"][i], geoms[key][i])
                per_unit[i][f"hamming_{key}"] = d
                deltas.append(d)
            entry["delta_sum"] = sum(deltas)
            entry["delta_mean"] = sum(deltas) / len(deltas)
            entry["theta_mean"] = position_error(
                centroids["original"], centroids[key]
            )
        result["maps"][key] = entry
    result["units"] = [per_unit[i] for i in keep]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("plot_json", nargs="+", type=Path)
    parser.add_argument(
        "--metrics",
        choices=("miaji", "2018", "all"),
        default="miaji",
        help="measures of Miaji et al. 2026 (default), of the 2018 SI, or both",
    )
    parser.add_argument(
        "--unit",
        choices=("region", "polygon"),
        default="region",
        help="for the 2018 measures: score every region (default) or polygon",
    )
    parser.add_argument(
        "--segment",
        type=float,
        default=0.01,
        help="resample rings for Frechet and Hausdorff (0: vertices only)",
    )
    parser.add_argument(
        "--rotation",
        choices=("local", "global", "map"),
        help="also compute each Miaji measure after a rotation: the nearest "
        "local minimum from the angle 0, the global minimum, or the region's "
        "own rotation under the map",
    )
    parser.add_argument("--json", type=Path, help="write per-unit results")
    args = parser.parse_args()

    results = [
        measures(p, args.unit, args.metrics, args.segment, args.rotation)
        for p in args.plot_json
    ]

    def short(r):
        return Path(r["file"]).name.replace("_conformality_plot.json", "")

    if args.metrics in ("miaji", "all"):
        print(
            "| Map | Regions | SymDiff T0 | SymDiff T | Fréchet T0 | Fréchet T "
            "| Hausdorff T0 | Hausdorff T |"
        )
        print("|---|---|---|---|---|---|---|---|")
        for r in results:
            m = r["miaji"]["maps"]
            print(
                f"| {short(r)} | {r['miaji']['n_regions']} | "
                f"{m['t0']['symdiff']:.4f} | {m['t']['symdiff']:.4f} | "
                f"{m['t0']['frechet']:.4f} | {m['t']['frechet']:.4f} | "
                f"{m['t0']['hausdorff']:.4f} | {m['t']['hausdorff']:.4f} |"
            )
            if r["n_rings_repaired"]:
                print(f"  ({r['n_rings_repaired']} invalid rings repaired)")
        if args.rotation:
            print(
                f"\nWith rotations permitted ({args.rotation} minimum; "
                "median |angle| in degrees):\n"
            )
            print(
                "| Map | SymDiff T0 | SymDiff T | Fréchet T0 | Fréchet T "
                "| Hausdorff T0 | Hausdorff T | SymDiff |angle| T0 | T |"
            )
            print("|---|---|---|---|---|---|---|---|---|")
            for r in results:
                m = r["miaji"]["maps"]
                print(
                    f"| {short(r)} | "
                    f"{m['t0']['symdiff_rot']:.4f} | {m['t']['symdiff_rot']:.4f} | "
                    f"{m['t0']['frechet_rot']:.4f} | {m['t']['frechet_rot']:.4f} | "
                    f"{m['t0']['hausdorff_rot']:.4f} | {m['t']['hausdorff_rot']:.4f} | "
                    f"{m['t0']['symdiff_median_abs_angle']:.1f} | "
                    f"{m['t']['symdiff_median_abs_angle']:.1f} |"
                )
    if args.metrics in ("2018", "all"):
        label = "Polygons" if args.unit == "polygon" else "Regions"
        print(f"| Map | {label} | α input | α T0 | α T | δ T0 | δ T | θ T0 | θ T |")
        print("|---|---|---|---|---|---|---|---|---|")
        for r in results:
            m = r["maps"]
            print(
                f"| {short(r)} | {r['n_units']} | {m['original']['alpha_mean']:.3f} | "
                f"{m['t0']['alpha_mean']:.3f} | {m['t']['alpha_mean']:.3f} | "
                f"{m['t0']['delta_sum']:.2f} | {m['t']['delta_sum']:.2f} | "
                f"{m['t0']['theta_mean']:.4f} | {m['t']['theta_mean']:.4f} |"
            )
            if r["n_rings_repaired"]:
                print(f"  ({r['n_rings_repaired']} invalid rings repaired)")
    if args.json:
        args.json.write_text(json.dumps(results, indent=1))
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
