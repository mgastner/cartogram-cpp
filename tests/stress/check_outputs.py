#!/usr/bin/env python3
"""Check the topology of the cartograms that the stress tests write.

Usage: check_outputs.py PATH [PATH ...]

Each PATH is a cartogram GeoJSON file or a folder, which is searched
recursively for files named *_cartogram.geojson. In each file:

- every polygon must be valid: its rings are simple, its holes lie inside
  its outer ring, and no hole lies inside another hole;
- no two polygons may overlap with positive area, whether they belong to
  different regions or to the same region.

Polygons may share edges. Neighboring regions share their borders, and the
parts of one region can share an edge, for example where the input splits a
region at the antimeridian and the projection joins the two halves again.

The exit code is 0 if every file passes, 1 if any file fails, and 2 if a
path does not exist or no file is found. Requires shapely 2.
"""

import json
import sys
from pathlib import Path

import shapely
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


def check(path):
    """Return the list of problems in one cartogram file"""
    with open(path, encoding="utf-8") as f:
        features = json.load(f)["features"]

    polygons = []
    regions = []
    for index, feature in enumerate(features):
        geometry = shape(feature["geometry"])
        for polygon in getattr(geometry, "geoms", [geometry]):
            polygons.append(polygon)
            regions.append(label(feature, index))

    problems = []
    for polygon, region in zip(polygons, regions):
        if not polygon.is_valid:
            problems.append(
                f"{region}: invalid polygon, "
                f"{shapely.is_valid_reason(polygon)}")

    total_area = sum(p.area for p in polygons)
    tree = shapely.STRtree(polygons)
    left, right = tree.query(polygons, predicate="intersects")
    for i, j in zip(left, right):
        if i >= j:
            continue
        area = shapely.intersection(polygons[i], polygons[j]).area
        if area > 0:
            pair = (
                f"two parts of {regions[i]}" if regions[i] == regions[j]
                else f"{regions[i]} and {regions[j]}")
            problems.append(
                f"{pair} overlap: area {area:.3g}, "
                f"{area / total_area:.3g} of the total area")
    return problems


def main():
    if len(sys.argv) < 2:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    files = cartogram_files(sys.argv[1:])
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
