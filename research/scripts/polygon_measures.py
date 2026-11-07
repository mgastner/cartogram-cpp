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
           [--segment LENGTH] [--json OUT.json]
"""

import argparse
import json
import math
import sys
from pathlib import Path

try:
    import numpy as np
    import shapely
    from shapely.affinity import scale, translate
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


def normalized(geom, region):
    """geom under the map that takes region to unit area with its centroid
    at the origin."""
    f = 1.0 / math.sqrt(region.area)
    c = region.centroid
    return scale(translate(geom, -c.x, -c.y), xfact=f, yfact=f, origin=(0, 0))


def miaji_measures(data, segment):
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
        units.append(entry)
    maps = {}
    for key in ("t0", "t"):
        maps[key] = {
            m: sum(u[f"{m}_{key}"] for u in units) / len(units)
            for m in ("symdiff", "frechet", "hausdorff", "hausdorff_sym")
        }
    return {"n_regions": len(units), "maps": maps, "units": units}, repaired


def measures(path, unit, metrics="2018", segment=0.01):
    data = json.loads(Path(path).read_text())
    if data["polygons"]["t"] is None:
        sys.exit(f"{path}: no corrected polygons (run with --conf_optimize)")
    result = {"file": str(path)}
    if metrics in ("miaji", "all"):
        result["miaji"], repaired = miaji_measures(data, segment)
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
    parser.add_argument("--json", type=Path, help="write per-unit results")
    args = parser.parse_args()

    results = [measures(p, args.unit, args.metrics, args.segment) for p in args.plot_json]

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
