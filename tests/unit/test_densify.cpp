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
// integration on a polygon with a thin part or a nearly degenerate crossing,
// with a uniform density, so that the map is the identity up to rounding,
// and the tests of needles also with nonuniform densities. The quadtree
// leaves are unit squares, so the triangulation has an edge on every lattice
// line and a diagonal in every lattice cell.

namespace
{

// Amplitudes A of the nonuniform densities
//   rho(i, j) = 1 + A cos(pi (i + 1/2) / l) cos(2 pi (j + 1/2) / l)
// on an l x l lattice, from the fourth review of PR #318. Each cell of the
// triangulation then has its own affine map, so that the projected image of
// a polygon edge bends where the edge crosses a triangle edge.
const std::vector<double> amplitudes = {1e-6, 1e-4, 1e-2, 0.1};

// Polygon with holes after one integration on an l x l lattice, with the
// density of amplitude A above, or with a uniform density if A is 0
Polygon_with_holes project_with_density(
  const Polygon_with_holes &in,
  const unsigned int l,
  const double amplitude = 0.0)
{
  Arguments args{};
  args.disable_simplification_densification = false;
  args.target_points_per_inset = 1000000;
  args.max_permitted_area_error = default_max_permitted_area_error;
  args.quadtree_leaf_count_factor = 1;

  InsetState inset("C", args);
  GeoDiv gd("polygon");
  gd.push_back(in);
  inset.push_back(gd);
  inset.set_grid_dimensions(l, l);
  FTReal2d &rho_init = inset.ref_to_rho_init();
  FTReal2d &rho_ft = inset.ref_to_rho_ft();
  rho_init.allocate(l, l);
  rho_ft.allocate(l, l);
  for (unsigned int i = 0; i < l; ++i) {
    for (unsigned int j = 0; j < l; ++j) {
      rho_init(i, j) = 1.0 + amplitude * std::cos(pi * (i + 0.5) / l) *
                               std::cos(2.0 * pi * (j + 0.5) / l);
      rho_ft(i, j) = 0.0;
    }
  }
  if (std::fpclassify(amplitude) == FP_ZERO) {
    rho_ft(0, 0) = 1.0;
  } else {
    inset.make_fftw_plans_for_rho();
    inset.execute_fftw_fwd_plan();
    inset.destroy_fftw_plans_for_rho();
  }
  inset.ref_to_fluxx_init().allocate(l, l);
  inset.ref_to_fluxy_init().allocate(l, l);
  inset.make_fftw_plans_for_flux();

  BOOST_TEST_REQUIRE(inset.flatten_density());
  BOOST_TEST_REQUIRE(inset.project());
  const Polygon_with_holes out = inset.geo_divs()[0].polygons_with_holes()[0];
  inset.destroy_fftw_plans_for_flux();
  inset.ref_to_fluxx_init().free();
  inset.ref_to_fluxy_init().free();
  rho_init.free();
  rho_ft.free();
  return out;
}

// Outer boundary of a polygon without holes after the same integration
Polygon project_with_density(
  const Polygon &in,
  const unsigned int l,
  const double amplitude = 0.0)
{
  return project_with_density(Polygon_with_holes(in), l, amplitude)
    .outer_boundary();
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

// A ring after the same integration, either as the outer boundary of a
// polygon without holes or as a hole in a square of side 8 around the center
// of the lattice, such as [60, 68] x [60, 68] for l = 128
Polygon project_ring_with_density(
  const Polygon &ring,
  const bool as_hole,
  const unsigned int l,
  const double amplitude = 0.0)
{
  if (!as_hole) {
    return project_with_density(ring, l, amplitude);
  }
  const double lo = 0.5 * l - 4.0;
  const double hi = 0.5 * l + 4.0;
  const Polygon outer = polygon_from({
    {lo, lo},
    {hi, lo},
    {hi, hi},
    {lo, hi},
  });
  Polygon hole = ring;
  if (!hole.is_clockwise_oriented()) {
    hole.reverse_orientation();
  }
  const std::vector<Polygon> holes = {hole};
  const Polygon_with_holes out = project_with_density(
    Polygon_with_holes(outer, holes.begin(), holes.end()),
    l,
    amplitude);
  BOOST_TEST_REQUIRE(out.number_of_holes() == 1u);
  BOOST_TEST(out.outer_boundary().is_simple());
  return *out.holes_begin();
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
  const Polygon out = project_with_density(in, 32);
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
// crossings, rounded to doubles, make the first needle cross itself. With the
// nonuniform densities, the first needle got the same point on both sides,
// and the projection put the second needle's tip across its other side.
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
    const Polygon out = project_with_density(in, l);
    BOOST_TEST(out.is_simple());
    BOOST_TEST(is_simple_with_exact_predicates(out));
    BOOST_TEST(largest_distance_from_boundary(out, in) < 1e-12);
    for (const double amplitude : amplitudes) {
      const Polygon out_nonuniform = project_with_density(in, l, amplitude);
      BOOST_TEST(out_nonuniform.is_simple());
      BOOST_TEST(is_simple_with_exact_predicates(out_nonuniform));
    }
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
  const Polygon out = project_with_density(in, 32);
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
// boundary is not simple". One of the two points is now moved by a unit in
// the last place.
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
  const Polygon out = project_with_density(in, 128);
  BOOST_TEST(out.is_simple());
  BOOST_TEST(is_simple_with_exact_predicates(out));
  BOOST_TEST(largest_distance_from_boundary(out, in) < 1e-12);
  for (const double amplitude : amplitudes) {
    const Polygon out_nonuniform = project_with_density(in, 128, amplitude);
    BOOST_TEST(out_nonuniform.is_simple());
    BOOST_TEST(is_simple_with_exact_predicates(out_nonuniform));
  }
}

// Holes from the second review of PR #318: needles like those above, whose
// two long sides are parallel and 5e-15 apart. Their ends lie 1.4e-8 from the
// diagonal x + y = 128, just farther than dbl_resolution, so both sides get a
// point where they cross the diagonal. The exact crossings lie 1.3e-7 apart
// along the needle. Rounding moves the crossing of one side 5e-15 off that
// side, which puts it 4e-22 beyond the other side, so that the hole crossed
// itself and the program stopped with "Hole is not simple". The crossing is
// now moved by a unit in the last place. Removing it instead, as an earlier
// revision of PR #318 did, made the second hole cross itself in the
// projection with the density of amplitude 1e-2.
BOOST_AUTO_TEST_CASE(holes_with_needles_along_a_diagonal_stay_simple)
{
  const std::vector<Polygon> holes = {
    polygon_from({
      {64.25000001, 63.75000001},
      {64.74999999, 63.24999999},
      {64.74999999, 63.24999998999999},
      {64.25000001, 63.750000009999994},
      {63.750000009999994, 63.25000001},
      {63.750000009999994, 64.25000001000001},
    }),
    polygon_from({
      {64.25000001, 63.75000001},
      {64.74999999, 63.24999999},
      {64.74999999, 63.249999990000006},
      {64.25000001, 63.75000001000001},
      {63.750000009999994, 64.25000001000001},
      {63.750000009999994, 63.25000001},
    }),
  };
  for (const Polygon &hole : holes) {
    BOOST_TEST_REQUIRE(is_simple_with_exact_predicates(hole));
    const Polygon out = project_ring_with_density(hole, true, 128);
    BOOST_TEST(out.is_simple());
    BOOST_TEST(is_simple_with_exact_predicates(out));
    BOOST_TEST(largest_distance_from_boundary(out, hole) < 1e-12);
    for (const double amplitude : amplitudes) {
      const Polygon out_nonuniform =
        project_ring_with_density(hole, true, 128, amplitude);
      BOOST_TEST(out_nonuniform.is_simple());
      BOOST_TEST(is_simple_with_exact_predicates(out_nonuniform));
    }
  }
}

// Needle from the third review of PR #318: its two sides are parallel to the
// diagonal x = y and a unit in the last place apart, and one side has a vertex
// at (64.5, 64.5). The other side crosses the diagonal x + y = 129 half a unit
// in the last place from (64.5, 64.5) in x and in y, and the crossing rounds
// to that point, so that the ring touched itself. The repair told inserted
// points from vertices by their coordinates, kept both occurrences, and the
// program stopped with "Outer boundary is not simple" or "Hole is not
// simple". Densification now records which points it inserted, and the
// inserted occurrence is moved by a unit in the last place.
//
// In the fourth review, a revision that removed the inserted occurrence
// failed with the nonuniform densities: without the point, the projection
// mapped the side as a straight segment between the images of its ends,
// whereas the image of the other side bends at (64.5, 64.5), so that the two
// sides crossed. With the density of amplitude 1e-4, the projection also
// rounded the image of the moved point onto the image of (64.5, 64.5).
BOOST_AUTO_TEST_CASE(point_inserted_at_a_vertex_of_the_other_side_is_moved)
{
  const Polygon in = polygon_from({
    {64.25, 64.25},
    {64.5, 64.5},
    {64.75, 64.75},
    {64.75000000000001, 64.75},
    {64.25000000000001, 64.25},
    {64.55, 63.95},
    {63.95, 63.95},
  });
  BOOST_TEST_REQUIRE(is_simple_with_exact_predicates(in));
  for (const bool as_hole : {false, true}) {
    const Polygon out = project_ring_with_density(in, as_hole, 128);
    BOOST_TEST(out.is_simple());
    BOOST_TEST(is_simple_with_exact_predicates(out));
    BOOST_TEST(largest_distance_from_boundary(out, in) < 1e-12);

    // The vertex stays where it is, and the other side keeps a point next to
    // it on the diagonal
    int at_vertex = 0;
    int near_vertex = 0;
    for (const Point &p : out) {
      if (p == Point(64.5, 64.5)) {
        ++at_vertex;
      }
      if (std::hypot(p.x() - 64.5, p.y() - 64.5) < 1e-12) {
        ++near_vertex;
      }
    }
    BOOST_TEST(at_vertex == 1);
    BOOST_TEST(near_vertex == 2);
    for (const double amplitude : amplitudes) {
      const Polygon out_nonuniform =
        project_ring_with_density(in, as_hole, 128, amplitude);
      BOOST_TEST(out_nonuniform.is_simple());
      BOOST_TEST(is_simple_with_exact_predicates(out_nonuniform));
    }
  }
}

// Needle from the fuzz test of the reply to the fourth review of PR #318
// (seed 2, case 24161), on a 32 x 32 lattice with the density of amplitude
// 0.1. As a hole, its projected image is simple with exact predicates, but
// the simplicity test of the Polygon type, whose predicates round, rejected
// it, and the program stopped with "Hole is not simple". The repair and
// InsetState::is_simple() now decide with exact predicates (fifth review).
BOOST_AUTO_TEST_CASE(hole_that_only_rounded_predicates_reject_is_accepted)
{
  const Polygon in = polygon_from({
    {16.530534662738543, 16.046951477438142},
    {16.741793058845442, 16.258206941154558},
    {16.853536754727838, 16.369949085968894},
    {16.853536754727838, 16.36994908596889},
    {16.530534662738543, 16.04695147743814},
    {16.530531718201015, 15.622687408736432},
    {16.106270594036832, 16.04695442197567},
  });
  BOOST_TEST_REQUIRE(is_simple_with_exact_predicates(in));
  for (const bool as_hole : {false, true}) {
    const Polygon out = project_ring_with_density(in, as_hole, 32, 0.1);
    BOOST_TEST(is_simple_with_exact_predicates(out));
  }
}
