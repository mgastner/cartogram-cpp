#define BOOST_TEST_MODULE test_topology_repair
#include "inset_state.hpp"
#include "topology_repair.hpp"
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <set>
#include <utility>
#include <vector>

namespace
{

const TopologyRepairSettings default_settings{1e-3, -1.0};

GeoDiv make_geodiv(
  const std::string &id,
  const std::vector<std::pair<double, double>> &ring)
{
  Polygon p;
  for (const auto &[x, y] : ring) {
    p.push_back(Point(x, y));
  }
  GeoDiv gd(id);
  gd.push_back(Polygon_with_holes(p));
  return gd;
}

std::vector<Point> vertices(const GeoDiv &gd)
{
  std::vector<Point> v;
  for (const auto &pwh : gd.polygons_with_holes()) {
    v.insert(
      v.end(),
      pwh.outer_boundary().vertices_begin(),
      pwh.outer_boundary().vertices_end());
    for (const auto &h : pwh.holes()) {
      v.insert(v.end(), h.vertices_begin(), h.vertices_end());
    }
  }
  return v;
}

double area(const GeoDiv &gd)
{
  double a = 0.0;
  for (const auto &pwh : gd.polygons_with_holes()) {
    a += pwh.outer_boundary().area();
    for (const auto &h : pwh.holes()) {
      a += h.area();
    }
  }
  return a;
}

// Vertices of the GeoDiv within distance tol of the segment from a to b
std::set<std::pair<double, double>> vertices_near(
  const GeoDiv &gd,
  const Point &a,
  const Point &b,
  const double tol)
{
  std::set<std::pair<double, double>> out;
  const Segment s(a, b);
  for (const auto &p : vertices(gd)) {
    if (CGAL::squared_distance(p, s) <= tol * tol) {
      out.insert({p.x(), p.y()});
    }
  }
  return out;
}

// Two triangles that share the border from (0, 0) to (10, 7). Polygon 2 has
// 99 extra vertices (0.1 k, 0.07 k) computed in double precision, most of
// which lie off the border by about 1e-16.
std::vector<GeoDiv> slanted_border()
{
  std::vector<std::pair<double, double>> p2 = {{0, 0}, {10, -3}, {10, 7}};
  for (int k = 99; k >= 1; --k) {
    p2.emplace_back(0.1 * k, 0.07 * k);
  }
  return {make_geodiv("A", {{0, 0}, {10, 7}, {0, 10}}), make_geodiv("B", p2)};
}

}  // namespace

BOOST_AUTO_TEST_CASE(clean_coverage_is_left_unchanged)
{
  std::vector<GeoDiv> gds = {
    make_geodiv("A", {{0, 0}, {10, 0}, {10, 5}, {10, 10}, {0, 10}}),
    make_geodiv("B", {{10, 0}, {20, 0}, {20, 10}, {10, 10}, {10, 5}})};
  const auto before_a = vertices(gds[0]);
  const auto before_b = vertices(gds[1]);
  const auto result = repair_topology(gds, default_settings);
  BOOST_TEST(result.within_tolerance);
  BOOST_TEST(result.changed.empty());
  BOOST_TEST((vertices(gds[0]) == before_a));
  BOOST_TEST((vertices(gds[1]) == before_b));
}

// The top edge has a small loop: the path crosses the edge at (5, 10)
BOOST_AUTO_TEST_CASE(small_self_intersection_is_repaired)
{
  std::vector<GeoDiv> gds = {make_geodiv(
    "A",
    {{0, 0}, {10, 0}, {10, 10}, {4, 10}, {5, 10.001}, {5, 9.999}, {0, 10}})};
  BOOST_TEST(!all_rings_simple(gds[0]));
  const auto result = repair_topology(gds, default_settings);
  BOOST_TEST(result.within_tolerance);
  BOOST_TEST_REQUIRE(result.changed.size() == 1u);
  BOOST_TEST(result.changed[0].had_self_intersection);
  BOOST_TEST(result.changed[0].relative_area_change < 1e-4);
  BOOST_TEST(all_rings_simple(gds[0]));
  BOOST_TEST(std::abs(area(gds[0]) - 100.0) < 1e-2);
}

BOOST_AUTO_TEST_CASE(repair_beyond_threshold_leaves_input_unchanged)
{
  std::vector<GeoDiv> gds = {make_geodiv(
    "A",
    {{0, 0}, {10, 0}, {10, 10}, {4, 10}, {5, 10.001}, {5, 9.999}, {0, 10}})};
  const auto before = vertices(gds[0]);
  const auto result = repair_topology(gds, {1e-9, -1.0});
  BOOST_TEST(!result.within_tolerance);
  BOOST_TEST((vertices(gds[0]) == before));
}

// A bow tie encloses two triangles whose signed areas cancel, so the repair
// changes the area completely
BOOST_AUTO_TEST_CASE(bow_tie_is_rejected)
{
  std::vector<GeoDiv> gds = {
    make_geodiv("A", {{0, 0}, {10, 10}, {10, 0}, {0, 10}})};
  const auto result = repair_topology(gds, default_settings);
  BOOST_TEST(!result.within_tolerance);
}

// B has a vertex 1e-12 off A's border; after the repair both regions have
// the same vertices along the border
BOOST_AUTO_TEST_CASE(near_miss_vertex_is_shared)
{
  std::vector<GeoDiv> gds = {
    make_geodiv("A", {{0, 0}, {10, 0}, {10, 10}, {0, 10}}),
    make_geodiv("B", {{10, 0}, {20, 0}, {20, 10}, {10, 10}, {10 + 1e-12, 5}})};
  const auto result = repair_topology(gds, default_settings);
  BOOST_TEST(result.within_tolerance);
  BOOST_TEST(!result.changed.empty());
  const auto a = vertices_near(gds[0], Point(10, 0), Point(10, 10), 1e-9);
  const auto b = vertices_near(gds[1], Point(10, 0), Point(10, 10), 1e-9);
  BOOST_TEST(a.size() == 3u);
  BOOST_TEST((a == b));
}

// B overlaps A in a strip of width 1e-7; the strip is given to one of them
BOOST_AUTO_TEST_CASE(overlap_is_merged)
{
  std::vector<GeoDiv> gds = {
    make_geodiv("A", {{0, 0}, {10, 0}, {10, 10}, {0, 10}}),
    make_geodiv("B", {{10 - 1e-7, 0}, {20, 0}, {20, 10}, {10 - 1e-7, 10}})};
  const auto result =
    repair_topology(gds, {1e-3, 1e-9});  // snapping alone cannot close it
  BOOST_TEST(result.within_tolerance);
  BOOST_TEST(std::abs(area(gds[0]) + area(gds[1]) - 200.0) < 1e-9);
}

// The shared slanted border keeps identical vertices in both regions after
// the repair and after simplification of the inset
BOOST_AUTO_TEST_CASE(slanted_shared_border_is_simplified_consistently)
{
  auto gds = slanted_border();
  const auto result = repair_topology(gds, default_settings);
  BOOST_TEST(result.within_tolerance);
  const Point a(0, 0);
  const Point b(10, 7);
  BOOST_TEST(
    (vertices_near(gds[0], a, b, 1e-9) == vertices_near(gds[1], a, b, 1e-9)));

  Arguments args{};
  InsetState inset("C", args);
  for (const auto &gd : gds) {
    inset.push_back(gd);
  }
  inset.simplify(1);
  const auto &out = inset.geo_divs();
  const auto shared_a = vertices_near(out[0], a, b, 1e-9);
  const auto shared_b = vertices_near(out[1], a, b, 1e-9);
  BOOST_TEST(shared_a.size() >= 2u);
  BOOST_TEST((shared_a == shared_b));
}
