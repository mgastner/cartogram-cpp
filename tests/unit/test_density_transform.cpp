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

constexpr unsigned int lx = 16;
constexpr unsigned int ly = 8;
constexpr double a = 0.5;
constexpr double b = 0.25;

// Density at the center of cell (i, j): a constant plus one cosine mode
// along each axis
double density(const unsigned int i, const unsigned int j)
{
  return 1.0 + a * std::cos(pi * (i + 0.5) / lx) +
         b * std::cos(3.0 * pi * (j + 0.5) / ly);
}

struct DensityTransform {
  Arguments args{};
  InsetState inset{"C", args};

  DensityTransform()
  {
    inset.set_grid_dimensions(lx, ly);
    inset.ref_to_rho_init().allocate(lx, ly);
    inset.ref_to_rho_ft().allocate(lx, ly);
    inset.make_fftw_plans_for_rho();
    for (unsigned int i = 0; i < lx; ++i) {
      for (unsigned int j = 0; j < ly; ++j) {
        inset.ref_to_rho_init()(i, j) = density(i, j);
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

}  // namespace

// In the normalization that flatten_density() expects, the inverse
// FFTW_REDFT01 transform weighs each nonzero frequency by 2. The coefficient
// of cos(pi m (i + 0.5) / lx) is therefore half its amplitude.
BOOST_FIXTURE_TEST_CASE(coefficients_of_known_density, DensityTransform)
{
  const FTReal2d &rho_ft = inset.ref_to_rho_ft();
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      double expected = 0.0;
      if (i == 0 && j == 0) {
        expected = 1.0;
      } else if (i == 1 && j == 0) {
        expected = a / 2.0;
      } else if (i == 0 && j == 3) {
        expected = b / 2.0;
      }
      BOOST_TEST(std::abs(rho_ft(i, j) - expected) < 1e-14);
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
      BOOST_TEST(std::abs(rho_init(i, j) - density(i, j)) < 1e-14);
    }
  }
}
