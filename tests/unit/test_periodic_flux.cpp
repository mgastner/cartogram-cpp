#define BOOST_TEST_MODULE test_periodic_flux
#include "constants.hpp"
#include "inset_state.hpp"
#include "interpolate_bilinearly.hpp"
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <random>
#include <vector>

// With --periodic_x, the density and the flux are periodic in x with period
// lx and have reflecting boundaries at y = 0 and y = ly. These tests check
// the flux that InsetState::compute_flux() derives from the density:
// - for a band-limited density with a known flux, the flux at the cell
//   centers agrees with the exact flux to rounding error;
// - the interpolated flux takes the same values at x and x + lx, and crosses
//   the seam at x = 0;
// - the y-component of the interpolated flux vanishes at y = 0 and y = ly;
// - for a random density, the spectral divergence of the flux equals
//   mean - rho to rounding error (except for the Nyquist frequency in x,
//   which the flux ignores), which is the continuity relation that the
//   velocity v = -F / rho(t) relies on;
// - the blur weights correspond to the wave numbers 2 pi k / lx.

namespace
{

// Term A cos(2 pi k x / lx + phase) cos(pi l y / ly) of the density
struct Mode {
  double amplitude;
  unsigned int k;
  double phase;
  unsigned int l;
};

const std::vector<Mode> modes = {
  {0.5, 1, 0.0, 0},
  {0.25, 3, -0.5 * pi, 2},
  {0.125, 0, 0.0, 3},
  {0.1, 2, 0.3, 1},
  {0.05, 5, 1.0, 7}};

Arguments periodic_arguments()
{
  Arguments args{};
  args.world = true;
  args.periodic_x = true;
  return args;
}

struct PeriodicInset {
  Arguments args = periodic_arguments();
  InsetState inset{"C", args};
  unsigned int lx;
  unsigned int ly;

  explicit PeriodicInset(
    const unsigned int grid_lx = 32,
    const unsigned int grid_ly = 16)
      : lx(grid_lx), ly(grid_ly)
  {
    inset.set_grid_dimensions(lx, ly);
    inset.ref_to_rho_init().allocate(lx, ly);
    inset.ref_to_rho_ft().allocate(lx, ly);
    inset.ref_to_fluxx_init().allocate(lx, ly);
    inset.ref_to_fluxy_init().allocate(lx, ly);
    inset.make_fftw_plans_for_rho();
    inset.make_fftw_plans_for_flux();
  }

  ~PeriodicInset()
  {
    inset.destroy_fftw_plans_for_rho();
    inset.destroy_fftw_plans_for_flux();
    inset.ref_to_rho_init().free();
    inset.ref_to_rho_ft().free();
    inset.ref_to_fluxx_init().free();
    inset.ref_to_fluxy_init().free();
  }

  double kx(const unsigned int k) const
  {
    return 2.0 * pi * k / lx;
  }

  double ky(const unsigned int l) const
  {
    return pi * l / ly;
  }

  double density(const double x, const double y) const
  {
    double rho = 1.0;
    for (const auto &m : modes) {
      rho +=
        m.amplitude * std::cos(kx(m.k) * x + m.phase) * std::cos(ky(m.l) * y);
    }
    return rho;
  }

  // Exact flux F = sum of A grad(mode) / K^2, which satisfies
  // div F = 1 - density
  double exact_flux(const double x, const double y, const char direction) const
  {
    double f = 0.0;
    for (const auto &m : modes) {
      const double k2 = kx(m.k) * kx(m.k) + ky(m.l) * ky(m.l);
      if (direction == 'x') {
        f += -m.amplitude * kx(m.k) * std::sin(kx(m.k) * x + m.phase) *
             std::cos(ky(m.l) * y) / k2;
      } else {
        f += -m.amplitude * ky(m.l) * std::cos(kx(m.k) * x + m.phase) *
             std::sin(ky(m.l) * y) / k2;
      }
    }
    return f;
  }

  void fill_with_manufactured_density()
  {
    for (unsigned int i = 0; i < lx; ++i) {
      for (unsigned int j = 0; j < ly; ++j) {
        inset.ref_to_rho_init()(i, j) = density(i + 0.5, j + 0.5);
      }
    }
  }

  // Flux component at cell (i, j), as flatten_density() reads it
  double flux(const unsigned int i, const unsigned int j, const char direction)
  {
    return (direction == 'x') ? inset.ref_to_fluxx_init()(i, j)
                              : inset.ref_to_fluxy_init()(i, j);
  }

  double interpolated_flux(const double x, const double y, const char d)
  {
    auto f = [&](unsigned int i, unsigned int j, char direction) {
      return flux(i, j, direction);
    };
    return interpolate_bilinearly_periodic_x(x, y, f, d, lx, ly);
  }
};

// blur_width() needs a grid whose long side has at least 256 cells
struct LargePeriodicInset : PeriodicInset {
  LargePeriodicInset() : PeriodicInset(256, 128) {}
};

}  // namespace

BOOST_FIXTURE_TEST_CASE(flux_of_band_limited_density_is_exact, PeriodicInset)
{
  fill_with_manufactured_density();
  inset.execute_fftw_fwd_plan();
  BOOST_TEST(std::abs(inset.ref_to_rho_ft()(0, 0) - 1.0) < 1e-14);
  inset.compute_flux();
  double max_err = 0.0;
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      for (const char d : {'x', 'y'}) {
        const double err =
          std::abs(flux(i, j, d) - exact_flux(i + 0.5, j + 0.5, d));
        max_err = std::max(max_err, err);
      }
    }
  }
  BOOST_TEST_MESSAGE("max |F - F_exact| = " << max_err);
  BOOST_TEST(max_err < 1e-13);
}

BOOST_FIXTURE_TEST_CASE(interpolated_flux_is_periodic_in_x, PeriodicInset)
{
  fill_with_manufactured_density();
  inset.execute_fftw_fwd_plan();
  inset.compute_flux();
  const double dlx = lx;
  for (const double x : {0.0, 0.2, 0.5, 3.7, 16.0, 31.5, 31.9}) {
    for (const double y : {0.0, 0.3, 4.5, 7.25, 15.8, 16.0}) {
      for (const char d : {'x', 'y'}) {
        const double f = interpolated_flux(x, y, d);
        BOOST_TEST_CONTEXT("x = " << x << ", y = " << y << ", d = " << d)
        {
          BOOST_TEST(std::abs(interpolated_flux(x + dlx, y, d) - f) < 1e-12);
          BOOST_TEST(std::abs(interpolated_flux(x - dlx, y, d) - f) < 1e-12);
          BOOST_TEST(
            std::abs(interpolated_flux(x + 3.0 * dlx, y, d) - f) < 1e-12);
        }
      }
    }
  }

  // At the seam, the interpolation averages the columns lx - 1 and 0, and
  // the x-component of the flux is nonzero: mass crosses the antimeridian
  for (unsigned int j = 0; j < ly; ++j) {
    const double y = j + 0.5;
    const double seam = interpolated_flux(0.0, y, 'x');
    const double average = 0.5 * (flux(lx - 1, j, 'x') + flux(0, j, 'x'));
    BOOST_TEST_CONTEXT("j = " << j)
    {
      BOOST_TEST(std::abs(seam - average) < 1e-14);
      BOOST_TEST(std::abs(interpolated_flux(dlx, y, 'x') - seam) < 1e-14);
    }
  }
  BOOST_TEST(std::abs(interpolated_flux(0.0, 0.5, 'x')) > 1e-2);
}

BOOST_FIXTURE_TEST_CASE(flux_satisfies_y_boundary_conditions, PeriodicInset)
{
  fill_with_manufactured_density();
  inset.execute_fftw_fwd_plan();
  inset.compute_flux();
  for (unsigned int n = 0; n <= 4 * lx; ++n) {
    const double x = 0.25 * n;
    BOOST_TEST_CONTEXT("x = " << x)
    {
      // No flux through the top and bottom sides
      BOOST_TEST(std::abs(interpolated_flux(x, 0.0, 'y')) < 1e-15);
      BOOST_TEST(
        std::abs(interpolated_flux(x, static_cast<double>(ly), 'y')) < 1e-15);

      // The exact flux has the same property, and the computed flux is
      // its sine series in y
      BOOST_TEST(std::abs(exact_flux(x, 0.0, 'y')) < 1e-15);
      BOOST_TEST(
        std::abs(exact_flux(x, static_cast<double>(ly), 'y')) < 1e-14);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(flux_satisfies_discrete_continuity, PeriodicInset)
{
  // Random density, which has content at every frequency
  std::mt19937 gen(20261115);
  std::uniform_real_distribution<double> dist(0.5, 1.5);
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      inset.ref_to_rho_init()(i, j) = dist(gen);
    }
  }
  inset.execute_fftw_fwd_plan();
  const FTReal2d &rho_ft = inset.ref_to_rho_ft();
  std::vector<double> rho_hat(lx * ly);
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      rho_hat[i * ly + j] = rho_ft(i, j);
    }
  }
  inset.compute_flux();

  // Normalized forward transforms of the flux components: F_x is a
  // halfcomplex series in x times a cosine series in y, and F_y is a
  // halfcomplex series in x times a sine series in y
  std::vector<double> fx_hat(lx * ly), fy_hat(lx * ly);
  fftw_plan px = fftw_plan_r2r_2d(
    static_cast<int>(lx),
    static_cast<int>(ly),
    inset.ref_to_fluxx_init().as_1d_array(),
    fx_hat.data(),
    FFTW_R2HC,
    FFTW_REDFT10,
    FFTW_ESTIMATE);
  fftw_plan py = fftw_plan_r2r_2d(
    static_cast<int>(lx),
    static_cast<int>(ly),
    inset.ref_to_fluxy_init().as_1d_array(),
    fy_hat.data(),
    FFTW_R2HC,
    FFTW_RODFT10,
    FFTW_ESTIMATE);
  fftw_execute(px);
  fftw_execute(py);
  fftw_destroy_plan(px);
  fftw_destroy_plan(py);
  const double norm = 1.0 / (2.0 * lx * ly);
  for (auto &c : fx_hat) {
    c *= norm;
  }
  for (auto &c : fy_hat) {
    c *= norm;
  }

  // Spectral divergence: d/dx multiplies frequency k by i k_x, and d/dy
  // turns sine coefficient j (frequency j + 1) into cosine coefficient j + 1
  // times k_y
  const unsigned int nyquist = lx / 2;
  double max_err = 0.0;
  double max_rho_hat = 0.0;
  for (unsigned int i = 0; i < lx; ++i) {
    const unsigned int k = std::min(i, lx - i);
    for (unsigned int j = 0; j < ly; ++j) {
      double div = 0.0;
      if (k != 0 && i != nyquist) {
        const double sign = (i < nyquist) ? -1.0 : 1.0;
        div += sign * kx(k) * fx_hat[(lx - i) * ly + j];
      }
      if (j > 0) {
        div += ky(j) * fy_hat[i * ly + (j - 1)];
      }

      // Expected: mean - rho, without the Nyquist frequency in x
      double expected = 0.0;
      if (i != nyquist && !(i == 0 && j == 0)) {
        expected = -rho_hat[i * ly + j];
      }
      max_err = std::max(max_err, std::abs(div - expected));
      max_rho_hat = std::max(max_rho_hat, std::abs(rho_hat[i * ly + j]));
    }

    // The sine series of F_y has no frequency ly
    max_err = std::max(max_err, std::abs(fy_hat[i * ly + (ly - 1)]));
  }
  BOOST_TEST_MESSAGE(
    "max |div F - (mean - rho)| in Fourier space = "
    << max_err << " (max |rho_hat| = " << max_rho_hat << ")");
  BOOST_TEST(max_rho_hat > 1e-2);
  BOOST_TEST(max_err < 1e-15);
}

// The blur must multiply frequency k of exp(2 pi i k x / lx) by the Gaussian
// weight of the wave number 2 pi k / lx
BOOST_FIXTURE_TEST_CASE(blur_uses_periodic_wave_numbers, LargePeriodicInset)
{
  const unsigned int k = 3;
  const double a = 0.5;
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      inset.ref_to_rho_init()(i, j) =
        1.0 + a * std::sin(kx(k) * (i + 0.5)) * std::cos(ky(2) * (j + 0.5));
    }
  }
  inset.execute_fftw_fwd_plan();
  const double bw = inset.blur_width();
  BOOST_TEST_REQUIRE(bw > 0.0);
  inset.blur_density();
  const double weight =
    std::exp(-0.5 * bw * bw * (kx(k) * kx(k) + ky(2) * ky(2)));
  double max_err = 0.0;
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      const double expected = 1.0 + a * weight * std::sin(kx(k) * (i + 0.5)) *
                                      std::cos(ky(2) * (j + 0.5));
      max_err =
        std::max(max_err, std::abs(inset.ref_to_rho_init()(i, j) - expected));
    }
  }
  BOOST_TEST_MESSAGE(
    "blur weight = " << weight << ", max error = " << max_err);
  BOOST_TEST(weight < 0.99);
  BOOST_TEST(max_err < 1e-14);
}
