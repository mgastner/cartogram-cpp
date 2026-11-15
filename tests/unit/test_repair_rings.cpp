#define BOOST_TEST_MODULE test_repair_rings
#include "repair_rings.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <vector>

// Rounding can make a ring cross or touch itself where another part of the
// ring is closer than a unit in the last place. The tests start from rings as
// rounding leaves them and check which points the repair moves or removes.

namespace
{

Polygon ring_from(const std::vector<std::pair<double, double>> &xy)
{
  Polygon p;
  for (const auto &[x, y] : xy) {
    p.push_back(Point(x, y));
  }
  return p;
}

// Simplicity with exact predicates, because the predicates of the
// project's Polygon type round
bool is_simple_with_exact_predicates(const Polygon &p)
{
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  CGAL::Polygon_2<K> q;
  for (const Point &v : p) {
    q.push_back(K::Point_2(v.x(), v.y()));
  }
  return q.is_simple();
}

bool contains(const Polygon &ring, const Point &p)
{
  return std::find(ring.begin(), ring.end(), p) != ring.end();
}

std::vector<GeoDiv> geo_divs_from(const std::vector<Polygon> &rings)
{
  std::vector<GeoDiv> geo_divs;
  for (const Polygon &ring : rings) {
    GeoDiv gd("G" + std::to_string(geo_divs.size()));
    gd.push_back(Polygon_with_holes(ring));
    geo_divs.push_back(gd);
  }
  return geo_divs;
}

const Polygon &ring_of(const std::vector<GeoDiv> &geo_divs, std::size_t g)
{
  return geo_divs[g].polygons_with_holes()[0].outer_boundary();
}

}  // namespace

// After the projection: the vertex Q lies one unit in the last place below
// the segment from (16, 16) to (20, 18), so that the segment from (18, 20) to
// Q crosses it. Q is also a vertex of a neighboring region, which shares the
// edge from Q to (16, 16). Only Q moves, in both rings, to the same point.
BOOST_AUTO_TEST_CASE(vertex_across_a_segment_moves_with_its_copy)
{
  const Point q(18.0, std::nextafter(17.0, 0.0));
  std::vector<GeoDiv> geo_divs = geo_divs_from({
    ring_from({{16.0, 16.0}, {20.0, 18.0}, {18.0, 20.0}, {q.x(), q.y()}}),
    ring_from({{16.0, 16.0}, {17.0, 12.0}, {19.0, 12.0}, {q.x(), q.y()}}),
  });
  BOOST_TEST_REQUIRE(!ring_of(geo_divs, 0).is_simple());
  BOOST_TEST_REQUIRE(ring_of(geo_divs, 1).is_simple());

  BOOST_TEST(repair_rings_after_projection(geo_divs));

  const Polygon &crossed = ring_of(geo_divs, 0);
  const Polygon &neighbor = ring_of(geo_divs, 1);
  BOOST_TEST(crossed.size() == 4u);
  BOOST_TEST(neighbor.size() == 4u);
  BOOST_TEST(crossed.is_simple());
  BOOST_TEST(is_simple_with_exact_predicates(crossed));
  BOOST_TEST(neighbor.is_simple());
  BOOST_TEST(!contains(crossed, q));
  BOOST_TEST(crossed[3] == neighbor[3]);
  BOOST_TEST(std::fabs(crossed[3].x() - q.x()) < 4e-15);
  BOOST_TEST(std::fabs(crossed[3].y() - q.y()) < 4e-15);
  for (const Point &p :
       {Point(16.0, 16.0), Point(20.0, 18.0), Point(18.0, 20.0)})
    BOOST_TEST(contains(crossed, p));
}

// After the projection: both sides of a needle one unit in the last place
// wide pass through the same point M, so that the ring touches itself there.
// One of the two occurrences of M moves, and the ring keeps its points.
BOOST_AUTO_TEST_CASE(point_where_a_ring_touches_itself_moves)
{
  const double u = std::nextafter(16.0, 17.0) - 16.0;
  std::vector<GeoDiv> geo_divs = geo_divs_from({ring_from({
    {16.0, 16.0},
    {16.5, 16.5},
    {17.0, 17.0},
    {17.0 + u, 17.0},
    {16.5, 16.5},
    {16.0 + u, 16.0},
  })});
  BOOST_TEST_REQUIRE(!ring_of(geo_divs, 0).is_simple());

  BOOST_TEST(repair_rings_after_projection(geo_divs));

  const Polygon &after = ring_of(geo_divs, 0);
  BOOST_TEST(after.size() == 6u);
  BOOST_TEST(after.is_simple());
  BOOST_TEST(is_simple_with_exact_predicates(after));
  BOOST_TEST(std::count(after.begin(), after.end(), Point(16.5, 16.5)) == 1);
}

// After densification: an inserted point lies far across the other side of
// the square, where no move by a unit in the last place helps. The point is
// removed, and the ring is the square from before densification again.
BOOST_AUTO_TEST_CASE(inserted_point_that_no_move_repairs_is_removed)
{
  const Polygon square =
    ring_from({{0.0, 0.0}, {4.0, 0.0}, {4.0, 4.0}, {0.0, 4.0}});
  const std::vector<GeoDiv> before = geo_divs_from({square});
  std::vector<GeoDiv> geo_divs = geo_divs_from({ring_from(
    {{0.0, 0.0}, {4.0, 0.0}, {-1.0, 2.0}, {4.0, 4.0}, {0.0, 4.0}})});
  std::vector<std::vector<bool>> inserted = {
    {false, false, true, false, false}};
  BOOST_TEST_REQUIRE(!ring_of(geo_divs, 0).is_simple());

  BOOST_TEST(repair_rings_after_densification(geo_divs, before, inserted));

  const Polygon &after = ring_of(geo_divs, 0);
  BOOST_TEST(after.size() == 4u);
  BOOST_TEST(std::equal(after.begin(), after.end(), square.begin()));
  BOOST_TEST(inserted[0].size() == 4u);
}
