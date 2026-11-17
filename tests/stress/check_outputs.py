#!/usr/bin/env python3
"""Check the topology of the cartograms that the stress tests write.

Usage: check_outputs.py PATH [PATH ...]

Each PATH is a cartogram GeoJSON file or a folder, which is searched
recursively for files named *_cartogram.geojson. Each file must be a
FeatureCollection with at least one feature, and:

- the geometry of every feature must be a nonempty Polygon or MultiPolygon
  whose rings are closed and have numeric coordinates;
- every polygon must be valid: its rings are simple, its holes lie inside
  its outer ring, and no hole lies inside another hole;
- no two valid polygons may overlap with positive area, whether they belong
  to different regions or to the same region.

Polygons may share edges. Neighboring regions share their borders, and the
parts of one region can share an edge, for example where the input splits a
region at the antimeridian and the projection joins the two halves again.

A problem in one file, including one that makes shapely raise an exception,
is reported for that file, and the other files are still checked. The exit
code is 0 if every file passes, 1 if any file fails, and 2 if no path is
given, a path does not exist, or no file is found. Requires shapely 2.
"""

import json
import numbers
import sys
from pathlib import Path

import shapely
from shapely.errors import GEOSException, ShapelyError
from shapely.geometry import shape


def cartogram_files(paths):
    """The files given and the cartogram files in the folders given, or None
    if a path does not exist"""
    files = []
    for path in map(Path, paths):
        if path.is_dir():
            files.extend(sorted(path.rglob("*_cartogram.geojson")))
        elif path.is_file():
            files.append(path)
        else:
            print(f"ERROR: {path} does not exist", file=sys.stderr)
            return None
    return files


def label(feature, index):
    """Region number and its first text property, for the report"""
    properties = feature.get("properties") or {}
    text = next((v for v in properties.values() if isinstance(v, str)), "")
    return f"region {index} ({text})" if text else f"region {index}"


def ring_problem(ring):
    """Why a ring of GeoJSON coordinates is malformed, or None"""
    if not isinstance(ring, list) or len(ring) < 4:
        return "a ring has fewer than 4 positions"
    for position in ring:
        if (not isinstance(position, list) or len(position) < 2
                or not all(isinstance(c, numbers.Real)
                           and not isinstance(c, bool)
                           for c in position[:2])):
            return "a ring has a position without two numeric coordinates"
    if ring[0][:2] != ring[-1][:2]:
        return "a ring is not closed"
    return None


def geometry_problem(geometry):
    """Why a feature's geometry is not a nonempty Polygon or MultiPolygon
    with well-formed rings, or None"""
    if not isinstance(geometry, dict):
        return "no geometry"
    kind = geometry.get("type")
    coordinates = geometry.get("coordinates")
    if kind == "Polygon":
        polygons = [coordinates]
    elif kind == "MultiPolygon":
        polygons = coordinates
    else:
        return f"geometry of type {kind}, not Polygon or MultiPolygon"
    if not isinstance(polygons, list) or not polygons:
        return f"empty {kind}"
    for rings in polygons:
        if not isinstance(rings, list) or not rings:
            return f"{kind} with an empty polygon"
        for ring in rings:
            problem = ring_problem(ring)
            if problem:
                return problem
    return None


def overlap_problems(polygons, regions):
    """Pairs of polygons that overlap with positive area"""
    problems = []
    if not polygons:
        return problems
    total_area = sum(p.area for p in polygons)
    tree = shapely.STRtree(polygons)
    left, right = tree.query(polygons, predicate="intersects")
    for i, j in zip(left, right):
        if i >= j:
            continue
        pair = (
            f"two parts of {regions[i]}" if regions[i] == regions[j]
            else f"{regions[i]} and {regions[j]}")
        try:
            area = shapely.intersection(polygons[i], polygons[j]).area
        except GEOSException as error:
            problems.append(f"{pair}: overlap test failed, {error}")
            continue
        if area > 0:
            problems.append(
                f"{pair} overlap: area {area:.3g}, "
                f"{area / total_area:.3g} of the total area")
    return problems


def check(path):
    """Return the list of problems in one cartogram file"""
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        return [f"cannot read the file: {error}"]
    features = data.get("features") if isinstance(data, dict) else None
    if not isinstance(features, list):
        return ["not a FeatureCollection"]
    if not features:
        return ["no features"]

    problems = []
    polygons = []
    regions = []
    for index, feature in enumerate(features):
        if not isinstance(feature, dict):
            problems.append(f"region {index}: not a feature")
            continue
        region = label(feature, index)
        problem = geometry_problem(feature.get("geometry"))
        if problem:
            problems.append(f"{region}: {problem}")
            continue
        try:
            geometry = shape(feature["geometry"])
            parts = list(getattr(geometry, "geoms", [geometry]))
        except (ValueError, TypeError, ShapelyError) as error:
            problems.append(f"{region}: cannot read the geometry, {error}")
            continue
        for polygon in parts:
            if polygon.is_empty:
                problems.append(f"{region}: empty polygon")
            elif not polygon.is_valid:
                problems.append(
                    f"{region}: invalid polygon, "
                    f"{shapely.is_valid_reason(polygon)}")
            else:
                polygons.append(polygon)
                regions.append(region)

    problems.extend(overlap_problems(polygons, regions))
    return problems


def main(arguments=None):
    arguments = sys.argv[1:] if arguments is None else arguments
    if not arguments:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    files = cartogram_files(arguments)
    if files is None:
        return 2
    if not files:
        print("ERROR: no cartogram GeoJSON files found", file=sys.stderr)
        return 2

    failed = 0
    for path in files:
        problems = check(path)
        status = "ok" if not problems else f"{len(problems)} problem(s)"
        print(f"{path}: {status}")
        for problem in problems:
            print(f"  {problem}")
        failed += bool(problems)

    print(f"{len(files) - failed} of {len(files)} files pass")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
