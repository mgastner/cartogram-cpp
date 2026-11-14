#define BOOST_TEST_MODULE test_densify
#include "constants.hpp"
#include "inset_state.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <limits>
#include <vector>

// Densification inserts a point wherever a polygon edge crosses an edge of
// the triangulation that projects the polygons. Each test runs one
// integration with a uniform density, so that the map is the identity up to
// rounding, on a polygon with a thin part or a nearly degenerate crossing.
// The quadtree leaves are unit squares, so the triangulation has an edge on
// every lattice line and a diagonal in every lattice cell.

namespace
{

// Outer boundary of the polygon after one integration with a uniform
// density on an l x l lattice
Polygon project_with_uniform_density(const Polygon &in, const unsigned int l)
{
  Arguments args{};
  args.disable_simplification_densification = false;
  args.target_points_per_inset = 1000000;
  args.max_permitted_area_error = default_max_permitted_area_error;
  args.quadtree_leaf_count_factor = 1;

  InsetState inset("C", args);
  GeoDiv gd("polygon");
  gd.push_back(Polygon_with_holes(in));
  inset.push_back(gd);
  inset.set_grid_dimensions(l, l);
  FTReal2d &rho_init = inset.ref_to_rho_init();
  FTReal2d &rho_ft = inset.ref_to_rho_ft();
  rho_init.allocate(l, l);
  rho_ft.allocate(l, l);
  for (unsigned int i = 0; i < l; ++i) {
    for (unsigned int j = 0; j < l; ++j) {
      rho_init(i, j) = 1.0;
      rho_ft(i, j) = 0.0;
    }
  }
  rho_ft(0, 0) = 1.0;
  inset.ref_to_fluxx_init().allocate(l, l);
  inset.ref_to_fluxy_init().allocate(l, l);
  inset.make_fftw_plans_for_flux();

  BOOST_TEST_REQUIRE(inset.flatten_density());
  BOOST_TEST_REQUIRE(inset.project());
  const Polygon out =
    inset.geo_divs()[0].polygons_with_holes()[0].outer_boundary();
  inset.destroy_fftw_plans_for_flux();
  inset.ref_to_fluxx_init().free();
  inset.ref_to_fluxy_init().free();
  rho_init.free();
  rho_ft.free();
  return out;
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

// Distance from p to the segment from a to b
double distance_to_segment(const Point &p, const Point &a, const Point &b)
{
  const double dx = b.x() - a.x();
  const double dy = b.y() - a.y();
  const double l2 = dx * dx + dy * dy;
  double t = ((p.x() - a.x()) * dx + (p.y() - a.y()) * dy) / l2;
  t = std::clamp(t, 0.0, 1.0);
  return std::hypot(p.x() - (a.x() + t * dx), p.y() - (a.y() + t * dy));
}

// Largest distance of a point of out from the boundary of in
double largest_distance_from_boundary(const Polygon &out, const Polygon &in)
{
  double largest = 0.0;
  for (const Point &p : out) {
    double d = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < in.size(); ++k) {
      d = std::min(d, distance_to_segment(p, in[k], in[(k + 1) % in.size()]));
    }
    largest = std::max(largest, d);
  }
  return largest;
}

// Polygon with the given vertices, in counterclockwise order
Polygon polygon_from(const std::vector<std::pair<double, double>> &xy)
{
  Polygon p;
  for (const auto &[x, y] : xy) {
    p.push_back(Point(x, y));
  }
  if (p.is_clockwise_oriented()) {
    p.reverse_orientation();
  }
  return p;
}

}  // namespace

// Thin band between y = 10.25 and y = 10.5. The top side crosses the
// lattice line x = 16 with an edge that is 9e-9 wide and 1e-9 high. Before
// PR #318, the edge counted as exactly vertical and received both ends of
// the vertical triangle edge on x = 16, which are a lattice unit apart; the
// resulting spike reached the bottom side of the band. The edge lies within
// dbl_resolution of x = 16 along its whole length, so it gets no point
// there.
BOOST_AUTO_TEST_CASE(short_edge_along_lattice_line_gets_no_point)
{
  const Polygon in = polygon_from({
    {4.0, 10.25},
    {28.0, 10.25},
    {28.0, 10.5},
    {16.0 + 5e-9, 10.5},
    {16.0 - 4e-9, 10.5 + 1e-9},
    {4.0, 10.5 + 1e-9},
  });
  const Polygon out = project_with_uniform_density(in, 32);
  BOOST_TEST(out.is_simple());
  BOOST_TEST(is_simple_with_exact_predicates(out));

  const double largest = largest_distance_from_boundary(out, in);
  BOOST_TEST_MESSAGE("largest distance from the band: " << largest);
  BOOST_TEST(largest < 1e-12);

  int on_line = 0;
  for (const Point &p : out) {
    if (!(p.x() < 16.0) && !(p.x() > 16.0) && p.y() > 10.4) {
      ++on_line;
    }
  }
  BOOST_TEST(on_line == 0);
}

// Needles from the review of PR #318: two sides one or two units in the
// last place apart, nearly parallel to the diagonal x + y = l of a lattice
// cell and about 1e-12 (or 1e-11) from it. Where the sides cross the
// diagonal is ill-conditioned: the two crossings lie about 3e-4 apart along
// the needle, and rounding the coordinates of an inserted point moves the
// point by about the width of the needle. The first commit of PR #318
// computed the crossing parameter with a relative error of 1e-3 and
// inserted the points, so that the needle crossed itself. Even the exact
// crossings, rounded to doubles, make the first needle cross itself.
BOOST_AUTO_TEST_CASE(needles_along_a_diagonal_stay_simple)
{
  const std::vector<std::pair<unsigned int, Polygon>> needles = {
    {32,
     polygon_from({
       {16.250000000001, 15.750000000001},
       {16.749999999999, 15.249999999999},
       {16.749999999999, 15.249999999999002},
       {16.250000000001, 15.750000000001002},
       {15.750000000000998, 16.250000000001002},
       {15.750000000000998, 15.250000000001},
     })},
    {32,
     polygon_from({
       {16.250000000001, 15.750000000001},
       {16.749999999999, 15.249999999999},
       {16.749999999999, 15.249999999999003},
       {16.250000000001, 15.750000000001004},
       {15.750000000000998, 16.250000000001002},
       {15.750000000000998, 15.250000000001},
     })},
    {128,
     polygon_from({
       {64.25000000001, 63.75000000001},
       {64.74999999999, 63.24999999999},
       {64.74999999999, 63.249999999989996},
       {64.25000000001, 63.75000000000999},
       {63.750000000010004, 63.25000000001},
       {63.750000000010004, 64.25000000001},
     })},
  };
  for (const auto &[l, in] : needles) {
    BOOST_TEST_REQUIRE(is_simple_with_exact_predicates(in));
    const Polygon out = project_with_uniform_density(in, l);
    BOOST_TEST(out.is_simple());
    BOOST_TEST(is_simple_with_exact_predicates(out));
    BOOST_TEST(largest_distance_from_boundary(out, in) < 1e-12);
  }
}

// An edge that crosses the diagonal x + y = 32 at a shallow angle, with its
// ends 2e-8 * sqrt(2) from the diagonal. The ends are symmetric about
// (16.5, 15.5) in floating point, so the crossing is exactly (16.5, 15.5).
// The first commit of PR #318 rounded ax + ay before subtracting it from 32
// and placed the crossing 1.6e-8 from (16.5, 15.5).
BOOST_AUTO_TEST_CASE(shallow_crossing_of_a_diagonal_is_accurate)
{
  const double ax = 16.25 + 2e-8;
  const double ay = 15.75 + 2e-8;
  const Polygon in = polygon_from({
    {ax, ay},
    {33.0 - ax, 31.0 - ay},
    {16.75, 15.75},
  });
  const Polygon out = project_with_uniform_density(in, 32);
  BOOST_TEST(out.is_simple());
  BOOST_TEST(is_simple_with_exact_predicates(out));

  double nearest = std::numeric_limits<double>::infinity();
  for (const Point &p : out) {
    nearest = std::min(nearest, std::hypot(p.x() - 16.5, p.y() - 15.5));
  }
  BOOST_TEST_MESSAGE(
    "distance of the nearest point from (16.5, 15.5): " << nearest);
  BOOST_TEST(nearest < 1e-13);
}

// Needle from a fuzz test: two sides one unit in the last place apart, at
// about 45 degrees, crossing the diagonals x + y = c at right angles. Where
// a diagonal crosses the needle, the exact crossings of the two sides are
// half a unit in the last place apart, and both round to the same point.
// The ring then touched itself, and the program stopped with "Outer
// boundary is not simple". Such points are now removed from the ring.
BOOST_AUTO_TEST_CASE(needle_across_a_diagonal_gets_no_shared_point)
{
  const Polygon in = polygon_from({
    {65.006382912452636, 63.634998335768685},
    {65.981577423198914, 64.610218542009576},
    {65.981577423198928, 64.610218542009576},
    {65.006382912452651, 63.634998335768685},
    {65.006388501866525, 63.210734267093578},
    {64.582118843777522, 63.634992746354804},
  });
  BOOST_TEST_REQUIRE(is_simple_with_exact_predicates(in));
  const Polygon out = project_with_uniform_density(in, 128);
  BOOST_TEST(out.is_simple());
  BOOST_TEST(is_simple_with_exact_predicates(out));
  BOOST_TEST(largest_distance_from_boundary(out, in) < 1e-12);
}
