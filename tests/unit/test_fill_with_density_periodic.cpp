#define BOOST_TEST_MODULE test_fill_with_density_periodic
#include "constants.hpp"
#include "inset_state.hpp"
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <vector>

// With --periodic_x, the density of a polygon that crosses the seam must be
// the density of the same polygon shifted into the interior of the lattice,
// shifted back by the same number of columns. The polygon below is a "C"
// that opens toward +x: right of the seam, its two prongs are disconnected,
// so the clipped copy that is rasterized has a bridge of zero width along
// the seam. No vertex lies on a grid line. With the vertices at integer
// coordinates, the rasterization of the polygon inside the lattice, which
// involves no seam, gave the ocean density to cells (15, 4) and (16, 4),
// although the polygon covers both cells completely.

namespace
{

constexpr unsigned int lx = 32;
constexpr unsigned int ly = 16;

GeoDiv c_shape(const double shift)
{
  Polygon p;
  for (const auto &[x, y] : std::vector<std::pair<double, double>>{
         {27.7, 4.3},
         {36.2, 4.3},
         {36.2, 6.4},
         {31.4, 6.4},
         {31.4, 9.6},
         {36.2, 9.6},
         {36.2, 11.7},
         {27.7, 11.7}}) {
    p.push_back(Point(x + shift, y));
  }
  GeoDiv gd("C");
  gd.push_back(Polygon_with_holes(p));
  return gd;
}

// Rasterized density of the C shifted by `shift` along x
std::vector<double> density(const double shift)
{
  Arguments args{};
  args.world = true;
  args.periodic_x = true;
  InsetState inset("C", args);
  const GeoDiv gd = c_shape(shift);
  inset.push_back(gd);
  inset.set_grid_dimensions(lx, ly);

  // A target area of twice the polygon area gives the polygon a nonzero
  // area error and, therefore, a nonzero weight in the density
  inset.insert_target_area("C", 2.0 * gd.area());
  inset.store_initial_area();
  inset.set_area_errors();
  inset.ref_to_rho_init().allocate(lx, ly);
  inset.ref_to_rho_ft().allocate(lx, ly);
  inset.make_fftw_plans_for_rho();
  inset.fill_with_density();
  std::vector<double> rho(lx * ly);
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      rho[i * ly + j] = inset.ref_to_rho_init()(i, j);
    }
  }

  // The polygons must be restored after the rasterization
  BOOST_TEST(std::abs(inset.geo_divs()[0].area() - gd.area()) < 1e-12);
  BOOST_TEST(
    inset.geo_divs()[0].polygons_with_holes()[0].outer_boundary().size() ==
    gd.polygons_with_holes()[0].outer_boundary().size());
  inset.destroy_fftw_plans_for_rho();
  inset.ref_to_rho_init().free();
  inset.ref_to_rho_ft().free();
  return rho;
}

}  // namespace

BOOST_AUTO_TEST_CASE(polygon_across_seam_has_shifted_density)
{
  // Across the seam: x from 27.7 to 36.2. Inside: x from 11.7 to 20.2.
  const std::vector<double> across = density(0.0);
  const std::vector<double> inside = density(-16.0);
  double max_diff = 0.0;
  double min_rho = dbl_inf;
  double max_rho = -dbl_inf;
  for (unsigned int i = 0; i < lx; ++i) {
    const unsigned int i_inside = (i + lx - 16) % lx;
    for (unsigned int j = 0; j < ly; ++j) {
      max_diff = std::max(
        max_diff,
        std::abs(across[i * ly + j] - inside[i_inside * ly + j]));
      min_rho = std::min(min_rho, across[i * ly + j]);
      max_rho = std::max(max_rho, across[i * ly + j]);
    }
  }
  BOOST_TEST_MESSAGE(
    "max difference = " << max_diff << ", density range = [" << min_rho << ", "
                        << max_rho << "]");

  // The polygon raises the density in the cells that it covers
  BOOST_TEST(max_rho - min_rho > 0.5);

  // The clipped copies lose a strip of width 1e-9 lx on each side of the
  // seam, so the densities differ slightly (measured: 1.1e-8)
  BOOST_TEST(max_diff < 1e-6);
}
