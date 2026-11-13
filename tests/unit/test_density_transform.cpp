#define BOOST_TEST_MODULE test_density_transform
#include "constants.hpp"
#include "inset_state.hpp"
#include <boost/test/included/unit_test.hpp>
#include <cmath>

// The forward transform of the density must be normalized whatever the blur
// width, because flatten_density() computes the flux from rho_ft_ and
// blur_density() leaves rho_ft_ unchanged when the blur width is 0.

namespace
{

constexpr double a = 0.5;
constexpr double b = 0.25;
constexpr double c = 0.125;

// Density at the center of cell (i, j) of an lx × ly grid: a constant, one
// cosine mode along each axis, and one product of modes along both axes
double density(
  const unsigned int i,
  const unsigned int j,
  const unsigned int lx,
  const unsigned int ly)
{
  const double x = (i + 0.5) / lx;
  const double y = (j + 0.5) / ly;
  return 1.0 + a * std::cos(pi * x) + b * std::cos(3.0 * pi * y) +
         c * std::cos(pi * x) * std::cos(2.0 * pi * y);
}

// Coefficient (i, j) of density() in the normalization that flatten_density()
// expects. The inverse FFTW_REDFT01 transform weights each nonzero frequency
// by 2 per axis, so the coefficient of a mode is its amplitude divided by 2
// for each axis along which the mode varies.
double expected_coefficient(const unsigned int i, const unsigned int j)
{
  if (i == 0 && j == 0) {
    return 1.0;
  }
  if (i == 1 && j == 0) {
    return a / 2.0;
  }
  if (i == 0 && j == 3) {
    return b / 2.0;
  }
  if (i == 1 && j == 2) {
    return c / 4.0;
  }
  return 0.0;
}

struct DensityTransform {
  Arguments args{};
  InsetState inset{"C", args};
  unsigned int lx;
  unsigned int ly;

  explicit DensityTransform(
    const unsigned int grid_lx = 16,
    const unsigned int grid_ly = 8)
      : lx(grid_lx), ly(grid_ly)
  {
    inset.set_grid_dimensions(lx, ly);
    inset.ref_to_rho_init().allocate(lx, ly);
    inset.ref_to_rho_ft().allocate(lx, ly);
    inset.make_fftw_plans_for_rho();
    for (unsigned int i = 0; i < lx; ++i) {
      for (unsigned int j = 0; j < ly; ++j) {
        inset.ref_to_rho_init()(i, j) = density(i, j, lx, ly);
      }
    }
    inset.execute_fftw_fwd_plan();
  }

  ~DensityTransform()
  {
    inset.destroy_fftw_plans_for_rho();
    inset.ref_to_rho_init().free();
    inset.ref_to_rho_ft().free();
  }
};

// blur_width() needs a grid whose long side has at least 256 cells
struct LargeDensityTransform : DensityTransform {
  LargeDensityTransform() : DensityTransform(256, 128) {}
};

}  // namespace

BOOST_FIXTURE_TEST_CASE(coefficients_of_known_density, DensityTransform)
{
  const FTReal2d &rho_ft = inset.ref_to_rho_ft();
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      BOOST_TEST_CONTEXT("i = " << i << ", j = " << j)
      {
        BOOST_TEST(
          std::abs(rho_ft(i, j) - expected_coefficient(i, j)) < 1e-14);
      }
    }
  }
}

// Without a blur, the backward transform must restore the density
BOOST_FIXTURE_TEST_CASE(backward_transform_restores_density, DensityTransform)
{
  inset.execute_fftw_bwd_plan();
  const FTReal2d &rho_init = inset.ref_to_rho_init();
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      BOOST_TEST_CONTEXT("i = " << i << ", j = " << j)
      {
        BOOST_TEST(std::abs(rho_init(i, j) - density(i, j, lx, ly)) < 1e-14);
      }
    }
  }
}

// The blur must apply only the Gaussian weights, not the normalization a
// second time: the mean density, rho_ft_(0, 0), stays unchanged, and every
// other coefficient shrinks by the weight of its frequency
BOOST_FIXTURE_TEST_CASE(
  blur_applies_only_the_gaussian_weights,
  LargeDensityTransform)
{
  const double bw = inset.blur_width();
  BOOST_TEST_REQUIRE(bw > 0.0);
  inset.blur_density();
  const FTReal2d &rho_ft = inset.ref_to_rho_ft();
  const double prefactor = -0.5 * bw * bw * pi * pi;
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      const double weight = std::exp(
        prefactor * (double(i) * i / (double(lx) * lx) +
                     double(j) * j / (double(ly) * ly)));
      BOOST_TEST_CONTEXT("i = " << i << ", j = " << j)
      {
        BOOST_TEST(
          std::abs(rho_ft(i, j) - weight * expected_coefficient(i, j)) <
          1e-14);
      }
    }
  }
}
