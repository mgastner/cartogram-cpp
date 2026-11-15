#define BOOST_TEST_MODULE test_remove_vertices_where_rings_cross
#include "inset_state.hpp"
#include "segments_that_meet.hpp"
#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <vector>

// Rounding in the projection can put a vertex across another segment of its
// ring when the vertex lies closer to the segment than a unit in the last
// place. The tests start from rings as rounding leaves them and check that
// InsetState::remove_vertices_where_rings_cross() removes the right points.

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

bool contains(const Polygon &ring, const Point &p)
{
  return std::find(ring.begin(), ring.end(), p) != ring.end();
}

}  // namespace

// The vertex Q lies one unit in the last place below the segment from A to
// B, so that the segment from C to Q crosses it. Q is also a vertex of a
// neighboring region. Only Q is removed, from both rings.
BOOST_AUTO_TEST_CASE(vertex_across_a_segment_is_removed_from_every_ring)
{
  const double q_y = std::nextafter(17.0, 0.0);
  const Polygon crossed =
    ring_from({{16.0, 16.0}, {20.0, 18.0}, {18.0, 20.0}, {18.0, q_y}});
  const Polygon neighbor =
    ring_from({{16.0, 16.0}, {17.0, 12.0}, {19.0, 12.0}, {18.0, q_y}});
  BOOST_TEST_REQUIRE(!crossed.is_simple());
  BOOST_TEST_REQUIRE(neighbor.is_simple());

  Arguments args{};
  InsetState inset("C", args);
  GeoDiv a("A");
  a.push_back(Polygon_with_holes(crossed));
  GeoDiv b("B");
  b.push_back(Polygon_with_holes(neighbor));
  inset.push_back(a);
  inset.push_back(b);

  BOOST_TEST(inset.remove_vertices_where_rings_cross());

  const Point q(18.0, q_y);
  const Polygon &crossed_after =
    inset.geo_divs()[0].polygons_with_holes()[0].outer_boundary();
  const Polygon &neighbor_after =
    inset.geo_divs()[1].polygons_with_holes()[0].outer_boundary();
  BOOST_TEST(crossed_after.size() == 3u);
  BOOST_TEST(neighbor_after.size() == 3u);
  BOOST_TEST(!contains(crossed_after, q));
  BOOST_TEST(!contains(neighbor_after, q));
  BOOST_TEST(crossed_after.is_simple());
  BOOST_TEST(is_simple_with_exact_predicates(crossed_after));
  BOOST_TEST(neighbor_after.is_simple());
}

// A spike of zero width: the ring runs up to (18, 23) and back down the same
// line to (18, 21). The segments next to the tip fold back, so the tip is
// removed, and (18, 21) lies on the segment up to the tip, so (18, 21) is
// removed as well. The spike is gone, and the other points remain.
BOOST_AUTO_TEST_CASE(tip_of_a_spike_of_zero_width_is_removed)
{
  const Polygon spike = ring_from(
    {{16.0, 16.0},
     {20.0, 16.0},
     {20.0, 20.0},
     {18.0, 20.0},
     {18.0, 23.0},
     {18.0, 21.0},
     {16.0, 20.0}});
  BOOST_TEST_REQUIRE(!spike.is_simple());

  Arguments args{};
  InsetState inset("C", args);
  GeoDiv a("A");
  a.push_back(Polygon_with_holes(spike));
  inset.push_back(a);

  BOOST_TEST(inset.remove_vertices_where_rings_cross());

  const Polygon &after =
    inset.geo_divs()[0].polygons_with_holes()[0].outer_boundary();
  BOOST_TEST(after.size() == spike.size() - 2);
  BOOST_TEST(!contains(after, Point(18.0, 23.0)));
  BOOST_TEST(!contains(after, Point(18.0, 21.0)));
  BOOST_TEST(after.is_simple());
}

// A degenerate triangle folds back at (18, 18). Removing that point would
// leave two points, so the ring keeps all three, and the function reports
// that a ring still fails.
BOOST_AUTO_TEST_CASE(ring_keeps_at_least_three_points)
{
  const Polygon degenerate =
    ring_from({{16.0, 16.0}, {18.0, 18.0}, {17.0, 17.0}});
  BOOST_TEST_REQUIRE(!degenerate.is_simple());

  Arguments args{};
  InsetState inset("C", args);
  GeoDiv a("A");
  a.push_back(Polygon_with_holes(degenerate));
  inset.push_back(a);

  BOOST_TEST(!inset.remove_vertices_where_rings_cross());
  BOOST_TEST(
    inset.geo_divs()[0].polygons_with_holes()[0].outer_boundary().size() ==
    3u);
}
