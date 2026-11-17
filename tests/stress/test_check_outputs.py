#!/usr/bin/env python3
"""Tests of check_outputs.py, run in CI before it checks the outputs of the
stress tests. Usage: test_check_outputs.py"""

import contextlib
import io
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

# Import the checker from this folder without leaving a __pycache__ folder
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_outputs  # noqa: E402


def rectangle(x0, y0, x1, y1):
    return [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]


def polygon(shell, *holes):
    return {"type": "Polygon", "coordinates": [shell, *holes]}


def multi(*polygons):
    return {"type": "MultiPolygon",
            "coordinates": [p["coordinates"] for p in polygons]}


def collection(*geometries):
    return {"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {"id": f"R{i}"}, "geometry": g}
        for i, g in enumerate(geometries)]}


SQUARE = polygon(rectangle(0, 0, 10, 10))
OVERLAP = polygon(rectangle(5, 0, 15, 10))
NEIGHBOR = polygon(rectangle(10, 0, 20, 10))
BOWTIE = polygon([[0, 0], [10, 10], [0, 10], [10, 0], [0, 0]])
WITH_HOLE = polygon(rectangle(0, 0, 20, 20), rectangle(5, 5, 15, 15))

# Name, contents of the file, expected exit code
CASES = [
    ("single_valid_polygon", collection(SQUARE), 0),
    ("shared_edge_between_regions", collection(SQUARE, NEIGHBOR), 0),
    ("shared_edge_within_region", collection(multi(SQUARE, NEIGHBOR)), 0),
    ("point_touch",
     collection(SQUARE, polygon(rectangle(10, 10, 20, 20))), 0),
    ("gap", collection(SQUARE, polygon(rectangle(10.01, 0, 20, 10))), 0),
    ("valid_hole", collection(WITH_HOLE), 0),
    ("island_in_hole",
     collection(WITH_HOLE, polygon(rectangle(7, 7, 13, 13))), 0),
    ("island_touching_hole_boundary",
     collection(WITH_HOLE, polygon(rectangle(5, 5, 15, 15))), 0),
    ("overlap_between_regions", collection(SQUARE, OVERLAP), 1),
    ("overlap_within_region", collection(multi(SQUARE, OVERLAP)), 1),
    ("complete_containment",
     collection(SQUARE, polygon(rectangle(2, 2, 8, 8))), 1),
    ("identical_polygons", collection(SQUARE, SQUARE), 1),
    ("one_ulp_overlap",
     collection(SQUARE,
                polygon(rectangle(math.nextafter(10, 0), 0, 20, 10))), 1),
    ("island_overlapping_hole_boundary",
     collection(WITH_HOLE, polygon(rectangle(4, 4, 15, 15))), 1),
    ("hole_outside_shell",
     collection(polygon(rectangle(0, 0, 10, 10),
                        rectangle(20, 20, 21, 21))), 1),
    ("nested_holes",
     collection(polygon(rectangle(0, 0, 20, 20), rectangle(2, 2, 18, 18),
                        rectangle(5, 5, 10, 10))), 1),
    ("self_intersection", collection(BOWTIE), 1),
    ("invalid_polygon_and_neighbor", collection(BOWTIE, SQUARE), 1),
    ("collapsed_polygon",
     collection(polygon([[0, 0], [1, 0], [2, 0], [0, 0]])), 1),
    ("null_coordinate",
     collection(polygon([[0, 0], [None, 0], [1, 1], [0, 0]])), 1),
    ("unclosed_ring",
     collection(polygon([[0, 0], [10, 0], [10, 10], [0, 10]])), 1),
    ("empty_feature_collection", collection(), 1),
    ("empty_polygon", collection({"type": "Polygon", "coordinates": []}), 1),
    ("empty_multipolygon",
     collection({"type": "MultiPolygon", "coordinates": []}), 1),
    ("point_only", collection({"type": "Point", "coordinates": [1, 1]}), 1),
    ("line_only",
     collection({"type": "LineString", "coordinates": [[0, 0], [1, 1]]}), 1),
    ("empty_geometry_collection",
     collection({"type": "GeometryCollection", "geometries": []}), 1),
    ("valid_polygon_and_empty_multipolygon",
     collection(SQUARE, {"type": "MultiPolygon", "coordinates": []}), 1),
    ("null_geometry", collection(None), 1),
]


def run(arguments):
    """Exit code and printed output of check_outputs.main()"""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = check_outputs.main([str(a) for a in arguments])
    return code, out.getvalue() + err.getvalue()


class CheckOutputsTest(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.root = Path(self.folder.name)

    def tearDown(self):
        self.folder.cleanup()

    def write(self, name, contents):
        path = self.root / f"{name}_cartogram.geojson"
        path.write_text(json.dumps(contents) if not isinstance(contents, str)
                        else contents)
        return path

    def test_cases(self):
        for name, contents, expected in CASES:
            with self.subTest(name):
                code, output = run([self.write(name, contents)])
                self.assertEqual(code, expected, output)
                self.assertNotIn("Traceback", output)

    def test_malformed_json(self):
        path = self.write("malformed", '{"type": "FeatureCollection", ')
        code, output = run([path])
        self.assertEqual(code, 1, output)
        self.assertIn("cannot read the file", output)

    def test_later_files_are_checked_after_a_failure(self):
        bad = self.write("a_invalid", collection(BOWTIE, SQUARE))
        good = self.write("b_valid", collection(SQUARE, NEIGHBOR))
        code, output = run([bad, good])
        self.assertEqual(code, 1, output)
        self.assertIn(f"{bad}: 1 problem(s)", output)
        self.assertIn(f"{good}: ok", output)
        self.assertIn("1 of 2 files pass", output)

    def test_folder_is_searched_recursively(self):
        sub = self.root / "sub"
        sub.mkdir()
        (sub / "x_cartogram.geojson").write_text(
            json.dumps(collection(SQUARE, OVERLAP)))
        self.write("y", collection(SQUARE))
        code, output = run([self.root])
        self.assertEqual(code, 1, output)
        self.assertIn("1 of 2 files pass", output)

    def test_usage_errors(self):
        empty = self.root / "empty"
        empty.mkdir()
        self.assertEqual(run([])[0], 2)
        self.assertEqual(run([self.root / "no_such_file"])[0], 2)
        self.assertEqual(run([empty])[0], 2)

    def test_file_and_empty_folder(self):
        empty = self.root / "empty"
        empty.mkdir()
        code, output = run([self.write("v", collection(SQUARE)), empty])
        self.assertEqual(code, 0, output)


if __name__ == "__main__":
    unittest.main()
