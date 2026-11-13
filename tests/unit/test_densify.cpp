#define BOOST_TEST_MODULE test_densify
#include "constants.hpp"
#include "inset_state.hpp"
#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <limits>

// Densification inserts a point wherever a polygon edge crosses an edge of
// the triangulation that projects the polygons. For an edge shorter than
// the tolerance dbl_resolution (1e-8 lattice units) in x, the old code
// treated the edge as exactly vertical. If it lay within the tolerance of a
// vertical triangle edge, the old code inserted both ends of that triangle
// edge, which are a lattice unit apart, so the polygon got a spike.
//
// This test runs one integration with a uniform density, so that the map is
// the identity, on a thin band whose top side crosses the lattice line
// x = 16 with such a short edge. The spike reached the bottom side of the
// band, and the program stopped with "Outer boundary is not simple".

namespace
{

constexpr unsigned int L = 32;

// Thin band between y = 10.25 and y = 10.5, in counterclockwise order. The
// top side crosses x = 16 with an edge that is 9e-9 wide and 1e-9 high.
Polygon band()
{
  Polygon p;
  p.push_back(Point(4.0, 10.25));
  p.push_back(Point(28.0, 10.25));
  p.push_back(Point(28.0, 10.5));
  p.push_back(Point(16.0 + 5e-9, 10.5));
  p.push_back(Point(16.0 - 4e-9, 10.5 + 1e-9));
  p.push_back(Point(4.0, 10.5 + 1e-9));
  return p;
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

// Run one integration with a uniform density and return the outer boundary
// of the projected band
Polygon project_band()
{
  Arguments args{};
  args.disable_simplification_densification = false;
  args.target_points_per_inset = 1000000;
  args.max_permitted_area_error = default_max_permitted_area_error;

  // Unit quadtree leaves: the triangulation has a vertical edge on every
  // lattice line x = i
  args.quadtree_leaf_count_factor = 1;

  InsetState inset("C", args);
  GeoDiv gd("band");
  gd.push_back(Polygon_with_holes(band()));
  inset.push_back(gd);
  inset.set_grid_dimensions(L, L);
  FTReal2d &rho_init = inset.ref_to_rho_init();
  FTReal2d &rho_ft = inset.ref_to_rho_ft();
  rho_init.allocate(L, L);
  rho_ft.allocate(L, L);
  for (unsigned int i = 0; i < L; ++i) {
    for (unsigned int j = 0; j < L; ++j) {
      rho_init(i, j) = 1.0;
      rho_ft(i, j) = 0.0;
    }
  }
  rho_ft(0, 0) = 1.0;
  inset.ref_to_fluxx_init().allocate(L, L);
  inset.ref_to_fluxy_init().allocate(L, L);
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

}  // namespace

BOOST_AUTO_TEST_CASE(short_edge_across_lattice_line_stays_on_its_segment)
{
  const Polygon in = band();
  const Polygon out = project_band();
  BOOST_TEST(out.is_simple());

  // Every point of the densified band lies on the original band (the map is
  // the identity up to rounding)
  double max_distance = 0.0;
  for (const Point &p : out) {
    double d = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < in.size(); ++k) {
      d = std::min(d, distance_to_segment(p, in[k], in[(k + 1) % in.size()]));
    }
    max_distance = std::max(max_distance, d);
  }
  BOOST_TEST_MESSAGE("largest distance from the band: " << max_distance);
  BOOST_TEST(max_distance < 1e-12);

  // The short edge gets exactly one point on x = 16
  int on_line = 0;
  for (const Point &p : out) {
    if (!(p.x() < 16.0) && !(p.x() > 16.0) && p.y() > 10.4) {
      ++on_line;
    }
  }
  BOOST_TEST(on_line == 1);
}
