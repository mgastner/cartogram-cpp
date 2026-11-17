#define BOOST_TEST_MODULE test_simplify_inset
#include "inset_state.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <boost/test/included/unit_test.hpp>
#include <set>
#include <utility>
#include <vector>

// InsetState::simplify() inserts all rings into one constrained triangulation.
// The tests check two cases of a border that two regions share: copies whose
// vertices lie exactly on each other's segments, and copies that only nearly
// coincide.

namespace
{

using XY = std::pair<double, double>;

// Points from p (inclusive) to q (exclusive) at a spacing of 1/4, so that each
// ring has more vertices than simplify() keeps
std::vector<XY> padding(const XY &p, const XY &q, const int n)
{
  std::vector<XY> points;
  for (int i = 0; i < n; ++i) {
    const double t = static_cast<double>(i) / n;
    points.emplace_back(
      p.first + t * (q.first - p.first),
      p.second + t * (q.second - p.second));
  }
  return points;
}

Polygon ring_from(const std::vector<XY> &xy)
{
  Polygon p;
  for (const auto &[x, y] : xy) {
    p.push_back(Point(x, y));
  }
  return p;
}

InsetState inset_from(const std::vector<std::vector<XY>> &rings)
{
  InsetState inset("C", Arguments{});
  for (const auto &xy : rings) {
    GeoDiv gd("G" + std::to_string(inset.n_geo_divs()));
    gd.push_back(Polygon_with_holes(ring_from(xy)));
    inset.push_back(gd);
  }
  return inset;
}

const Polygon &ring_of(const InsetState &inset, const std::size_t g)
{
  return inset.geo_divs()[g].polygons_with_holes()[0].outer_boundary();
}

bool simple_with_exact_predicates(const Polygon &p)
{
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  CGAL::Polygon_2<K> q;
  for (const Point &v : p) {
    q.push_back(K::Point_2(v.x(), v.y()));
  }
  return q.is_simple();
}

// Vertices of a ring on the vertical line x = 10
std::set<XY> on_line_x10(const Polygon &ring)
{
  std::set<XY> points;
  for (const Point &v : ring) {
    if (v.x() == 10.0) {
      points.emplace(v.x(), v.y());
    }
  }
  return points;
}

}  // namespace

// The square [0, 10] x [0, 10] and the rectangle [10, 19] x [0, 10] share the
// border x = 10. The two copies of the border have different vertices, all
// exactly on the line, so the triangulation merges the copies and both
// regions keep the same vertices on the border.
BOOST_AUTO_TEST_CASE(
  exactly_shared_border_with_different_vertices_stays_shared)
{
  std::vector<XY> left = padding({0, 0}, {10, 0}, 40);
  for (const double y : {0.0, 2.5, 5.0, 7.5}) {
    left.emplace_back(10, y);
  }
  for (const XY &p : padding({10, 10}, {0, 10}, 40)) {
    left.push_back(p);
  }
  for (const XY &p : padding({0, 10}, {0, 0}, 40)) {
    left.push_back(p);
  }
  std::vector<XY> right = padding({10, 0}, {19, 0}, 36);
  for (const XY &p : padding({19, 0}, {19, 10}, 40)) {
    right.push_back(p);
  }
  for (const XY &p : padding({19, 10}, {10, 10}, 36)) {
    right.push_back(p);
  }
  for (const double y : {10.0, 8.0, 6.0, 4.0, 2.0}) {
    right.emplace_back(10, y);
  }

  InsetState inset = inset_from({left, right});
  BOOST_TEST(
    (on_line_x10(ring_of(inset, 0)) != on_line_x10(ring_of(inset, 1))));
  inset.simplify(30);

  BOOST_TEST(
    (on_line_x10(ring_of(inset, 0)) == on_line_x10(ring_of(inset, 1))));
  BOOST_TEST(simple_with_exact_predicates(ring_of(inset, 0)));
  BOOST_TEST(simple_with_exact_predicates(ring_of(inset, 1)));
}

// The two regions share the border from (0, 0) to (10, 7). The copy of the
// second region has the points (k, 0.7 * k), which lie just below the border
// for k = 1, ..., 9, except k = 5. With plain floating-point predicates, the
// triangulation calls some of them collinear with the border, merges the
// copies in part, and constructs an intersection point, (5.5, 3.85), that
// neither input has. With exact predicates, simplify() keeps only input
// vertices.
BOOST_AUTO_TEST_CASE(nearly_collinear_border_keeps_only_input_vertices)
{
  std::vector<XY> upper = {{0, 0}, {10, 7}};
  for (const XY &p : padding({0, 10}, {0, 0}, 40)) {
    upper.push_back(p);
  }
  std::vector<XY> lower = {{0, 0}};
  for (const XY &p : padding({10, -3}, {10, 7}, 40)) {
    lower.push_back(p);
  }
  lower.emplace_back(10, 7);
  for (int k = 9; k >= 1; --k) {
    lower.emplace_back(k, 0.7 * k);
  }

  InsetState inset = inset_from({upper, lower});
  std::set<XY> input;
  for (std::size_t g = 0; g < 2; ++g) {
    for (const Point &v : ring_of(inset, g)) {
      input.emplace(v.x(), v.y());
    }
  }
  inset.simplify(30);

  for (std::size_t g = 0; g < 2; ++g) {
    for (const Point &v : ring_of(inset, g)) {
      BOOST_TEST(
        input.contains({v.x(), v.y()}),
        "vertex (" << v.x() << ", " << v.y() << ") of ring " << g
                   << " is not an input vertex");
    }
    BOOST_TEST(simple_with_exact_predicates(ring_of(inset, g)));
  }
}
