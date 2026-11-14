#define BOOST_TEST_MODULE test_flatten_density_periodic
#include "constants.hpp"
#include "inset_state.hpp"
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <vector>

// End-to-end check of the flow-based method with --periodic_x on a density
// with a known answer. The density
//   rho(x, y) = 1 + a sin(2 pi x / L)
// depends on x only, so the flow moves points only in x. Its flux at the seam
// x = 0 is nonzero, so points cross the seam during the integration. Every
// point keeps the mass between itself and the seam, and the integrator
// translates the map so that the seam returns to x = 0. When the density is
// equalized at t = 1, a point that started at x has reached
//   X(x) = x + (a L / (2 pi)) (1 - cos(2 pi x / L)),
// which is the mass between 0 and x. We project the vertices of a polygon
// that crosses the seam and compare them with this map.

namespace
{

constexpr unsigned int L = 64;
constexpr double a = 0.5;

double exact_x(const double x)
{
  return x + a * L / (2.0 * pi) * (1.0 - std::cos(2.0 * pi * x / L));
}

// Rectangle [x0, x0 + L] x [y0, y1] with a vertex at every integer x on its
// horizontal sides, in counterclockwise order
GeoDiv strip(const int x0, const double y0, const double y1)
{
  Polygon p;
  for (int i = 0; i <= static_cast<int>(L); ++i) {
    p.push_back(Point(x0 + i, y0));
  }
  for (int i = static_cast<int>(L) + 1; i-- > 0;) {
    p.push_back(Point(x0 + i, y1));
  }
  GeoDiv gd("strip");
  gd.push_back(Polygon_with_holes(p));
  return gd;
}

// Run one integration with the flow-based method and return the projected
// vertices of the strip
std::vector<Point> project_strip(const GeoDiv &gd)
{
  Arguments args{};
  args.world = true;
  args.periodic_x = true;
  args.disable_simplification_densification = true;
  args.max_permitted_area_error = default_max_permitted_area_error;

  // Refine the quadtree down to unit cells so that every lattice point is a
  // quadtree corner and is moved by the integrator itself
  args.quadtree_leaf_count_factor = 1;

  InsetState inset("C", args);
  inset.push_back(gd);
  inset.set_grid_dimensions(L, L);
  FTReal2d &rho_init = inset.ref_to_rho_init();
  FTReal2d &rho_ft = inset.ref_to_rho_ft();
  rho_init.allocate(L, L);
  rho_ft.allocate(L, L);
  inset.make_fftw_plans_for_rho();
  for (unsigned int i = 0; i < L; ++i) {
    for (unsigned int j = 0; j < L; ++j) {
      rho_init(i, j) = 1.0 + a * std::sin(2.0 * pi * (i + 0.5) / L);
    }
  }
  inset.execute_fftw_fwd_plan();

  // The flux is computed with inverse Fourier transforms
  inset.ref_to_fluxx_init().allocate(L, L);
  inset.ref_to_fluxy_init().allocate(L, L);
  inset.make_fftw_plans_for_flux();

  BOOST_TEST_REQUIRE(inset.flatten_density());
  BOOST_TEST_REQUIRE(inset.project());
  const Polygon &out =
    inset.geo_divs()[0].polygons_with_holes()[0].outer_boundary();
  std::vector<Point> v(out.vertices_begin(), out.vertices_end());
  inset.destroy_fftw_plans_for_flux();
  inset.destroy_fftw_plans_for_rho();
  inset.ref_to_fluxx_init().free();
  inset.ref_to_fluxy_init().free();
  rho_init.free();
  rho_ft.free();
  return v;
}

}  // namespace

// The strip covers [-16, 48], so that the projection must also locate points
// with x < 0 through their copies at x + L
BOOST_AUTO_TEST_CASE(sine_density_reaches_exact_map_across_the_seam)
{
  const GeoDiv gd = strip(-16, 8.0, 24.0);
  const Polygon &in = gd.polygons_with_holes()[0].outer_boundary();
  const std::vector<Point> out = project_strip(gd);
  BOOST_TEST_REQUIRE(out.size() == in.size());
  double max_dx = 0.0;
  double max_dy = 0.0;
  for (std::size_t k = 0; k < out.size(); ++k) {
    const Point &p = in[k];
    max_dx = std::max(max_dx, std::abs(out[k].x() - exact_x(p.x())));
    max_dy = std::max(max_dy, std::abs(out[k].y() - p.y()));
  }
  BOOST_TEST_MESSAGE("max |X - X_exact| = " << max_dx);

  // During the integration, the points move by up to a L / (2 pi) = 5.1
  // lattice units. The remaining error comes from the bilinear interpolation
  // of the velocity and from the time steps (measured: 8.1e-3 lattice units,
  // against 2.6e-3 in test_flatten_density, whose velocity has half the wave
  // number and twice the amplitude).
  BOOST_TEST(max_dx < 2e-2);

  // The velocity has no y component
  BOOST_TEST(max_dy < 1e-12);
}
