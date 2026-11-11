#define BOOST_TEST_MODULE test_flatten_density
#include "constants.hpp"
#include "inset_state.hpp"
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <vector>

// End-to-end check of the flow-based method on a density with a known
// answer. The density
//   rho(x, y) = 1 + a cos(pi x / L)
// depends on x only, so the flow moves points only in x. At time t, the
// interpolated density is rho(x, t) = 1 + (1 - t) a cos(pi x / L), and every
// point keeps the mass to its left. Hence a point that starts at x is at
// X(x, t), where
//   X + (1 - t) (a L / pi) sin(pi X / L) = x + (a L / pi) sin(pi x / L),
// and Y = y. This holds also for t > 1 (over-relaxation), as long as
// rho(x, t) > 0. At t = 1, when the density is equalized,
//   X(x, 1) = x + (a L / pi) sin(pi x / L).
// We project the vertices of a polygon that lie on the integer lattice and
// compare them with this map.

namespace
{

constexpr unsigned int L = 64;
constexpr double a = 0.5;

double exact_x(const double x, const double t = 1.0)
{
  const double rhs = x + a * L / pi * std::sin(pi * x / L);
  const double c = (1.0 - t) * a * L / pi;

  // Newton's method; the derivative 1 + (1 - t) a cos(pi X / L) is the
  // positive density rho(X, t)
  double xx = rhs;
  for (int k = 0; k < 50; ++k) {
    const double f = xx + c * std::sin(pi * xx / L) - rhs;
    xx -= f / (1.0 + c * pi / L * std::cos(pi * xx / L));
  }
  return xx;
}

// Rectangle [0, L] x [y0, y1] with a vertex at every integer x on its
// horizontal sides, in counterclockwise order
GeoDiv strip(const double y0, const double y1)
{
  Polygon p;
  for (unsigned int i = 0; i <= L; ++i) {
    p.push_back(Point(i, y0));
  }
  for (unsigned int i = L + 1; i-- > 0;) {
    p.push_back(Point(i, y1));
  }
  GeoDiv gd("strip");
  gd.push_back(Polygon_with_holes(p));
  return gd;
}

// Run one integration with the flow-based method, with over-relaxation
// omega, and return the projected vertices of the strip and the time reached
std::vector<Point> project_strip(
  const GeoDiv &gd,
  const double omega,
  double &t_reached)
{
  Arguments args{};
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
  for (unsigned int i = 0; i < L; ++i) {
    for (unsigned int j = 0; j < L; ++j) {
      rho_init(i, j) = 1.0 + a * std::cos(pi * (i + 0.5) / L);
      rho_ft(i, j) = 0.0;
    }
  }

  // Cosine coefficients in the normalization that flatten_density() expects:
  // the inverse FFTW_REDFT01 transform weighs frequency 1 by 2
  rho_ft(0, 0) = 1.0;
  rho_ft(1, 0) = a / 2.0;

  // The flux is computed with inverse Fourier transforms
  inset.ref_to_fluxx_init().allocate(L, L);
  inset.ref_to_fluxy_init().allocate(L, L);
  inset.make_fftw_plans_for_flux();

  inset.set_flow_overrelaxation(omega);
  BOOST_TEST_REQUIRE(inset.flatten_density());
  t_reached = inset.flow_t_reached();
  BOOST_TEST_REQUIRE(inset.project());
  const Polygon &out =
    inset.geo_divs()[0].polygons_with_holes()[0].outer_boundary();
  std::vector<Point> v(out.vertices_begin(), out.vertices_end());
  inset.destroy_fftw_plans_for_flux();
  inset.ref_to_fluxx_init().free();
  inset.ref_to_fluxy_init().free();
  rho_init.free();
  rho_ft.free();
  return v;
}

// Largest deviations of the projected strip from the exact map at time t
void compare_with_exact_map(
  const GeoDiv &gd,
  const std::vector<Point> &out,
  const double t,
  double &max_dx,
  double &max_dy)
{
  const Polygon &in = gd.polygons_with_holes()[0].outer_boundary();
  BOOST_TEST_REQUIRE(out.size() == in.size());
  max_dx = 0.0;
  max_dy = 0.0;
  for (std::size_t k = 0; k < out.size(); ++k) {
    const Point &p = in[k];
    max_dx = std::max(max_dx, std::abs(out[k].x() - exact_x(p.x(), t)));
    max_dy = std::max(max_dy, std::abs(out[k].y() - p.y()));
  }
}

}  // namespace

BOOST_AUTO_TEST_CASE(one_dimensional_density_reaches_exact_map)
{
  const GeoDiv gd = strip(8.0, 24.0);
  double t_reached = 0.0;
  const std::vector<Point> out = project_strip(gd, 0.0, t_reached);
  BOOST_TEST(t_reached == 1.0);
  double max_dx = 0.0;
  double max_dy = 0.0;
  compare_with_exact_map(gd, out, 1.0, max_dx, max_dy);
  BOOST_TEST_MESSAGE("max |X - X_exact| = " << max_dx);

  // The points move up to a L / pi = 10 lattice units. The remaining error
  // comes from the bilinear interpolation of the velocity and from the time
  // steps (measured: 2.6e-3 lattice units). When the last time step ended
  // past t = 1, the error was 0.53 lattice units.
  BOOST_TEST(max_dx < 1e-2);

  // The velocity has no y component
  BOOST_TEST(max_dy == 0.0);
}

BOOST_AUTO_TEST_CASE(overrelaxed_integration_reaches_exact_map)
{
  // With omega = 0.5, the integration ends at t = 1.5, where the interpolated
  // density 1 - 0.5 a cos(pi x / L) is at least 0.75
  const GeoDiv gd = strip(8.0, 24.0);
  double t_reached = 0.0;
  const std::vector<Point> out = project_strip(gd, 0.5, t_reached);
  BOOST_TEST(t_reached == 1.5);
  double max_dx = 0.0;
  double max_dy = 0.0;
  compare_with_exact_map(gd, out, 1.5, max_dx, max_dy);
  BOOST_TEST_MESSAGE("max |X - X_exact| at t = 1.5: " << max_dx);

  // Measured: 2.7e-3 lattice units
  BOOST_TEST(max_dx < 1e-2);
  BOOST_TEST(max_dy == 0.0);
}

BOOST_AUTO_TEST_CASE(overrelaxation_is_limited_by_density_maximum)
{
  // The largest grid value of rho_init is 1 + a cos(pi / (2 L)). The
  // interpolated density must stay at least min_flow_density_fraction of the
  // mean 1, so omega is at most
  //   (1 - min_flow_density_fraction) / (a cos(pi / (2 L))).
  const double omega_limit =
    (1.0 - min_flow_density_fraction) / (a * std::cos(pi / (2.0 * L)));
  const GeoDiv gd = strip(8.0, 24.0);
  double t_reached = 0.0;
  const std::vector<Point> out = project_strip(gd, 5.0, t_reached);
  BOOST_TEST(
    t_reached == 1.0 + omega_limit,
    boost::test_tools::tolerance(1e-12));
  double max_dx = 0.0;
  double max_dy = 0.0;
  compare_with_exact_map(gd, out, t_reached, max_dx, max_dy);
  BOOST_TEST_MESSAGE(
    "max |X - X_exact| at t = " << t_reached << ": " << max_dx);

  // Measured: 5.1e-3 lattice units at t = 2.0003
  BOOST_TEST(max_dx < 1e-2);
  BOOST_TEST(max_dy == 0.0);
}
