#define BOOST_TEST_MODULE test_conformality_spectral_basis
#include "conformality/spectral_basis.hpp"
#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <random>
#include <stdexcept>
#include <vector>

using conformality::fill_harmonics;
using conformality::Mat2;
using conformality::SpectralBasis;
using conformality::Vec2;

namespace
{

std::vector<double> random_coefficients(std::size_t n, std::mt19937 &rng)
{
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  std::vector<double> c(n);
  for (auto &x : c) {
    x = dist(rng);
  }
  return c;
}

// Upper bound on |c_mn| * kmax^p summed over modes, used to scale tolerances
double coefficient_scale(
  const SpectralBasis &basis,
  const std::vector<double> &c,
  int power)
{
  const double kmax = std::max(basis.kx(basis.mx()), basis.ky(basis.my()));
  double sum = 0.0;
  for (const double x : c) {
    sum += std::abs(x);
  }
  return sum * std::pow(kmax, power);
}

double frobenius_inner(const Mat2 &a, const Mat2 &b)
{
  return a.a11 * b.a11 + a.a12 * b.a12 + a.a21 * b.a21 + a.a22 * b.a22;
}

}  // namespace

BOOST_AUTO_TEST_CASE(constructor_rejects_invalid_arguments)
{
  BOOST_CHECK_THROW(SpectralBasis(0.0, 1.0, 2, 2), std::invalid_argument);
  BOOST_CHECK_THROW(SpectralBasis(1.0, -1.0, 2, 2), std::invalid_argument);
  BOOST_CHECK_THROW(SpectralBasis(1.0, 1.0, 0, 2), std::invalid_argument);
  BOOST_CHECK_THROW(SpectralBasis(1.0, 1.0, 2, 0), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(harmonics_match_standard_library)
{
  constexpr std::size_t count = 256;
  std::vector<double> s(count), c(count);
  std::mt19937 rng(12345);
  std::uniform_real_distribution<double> dist(0.0, std::numbers::pi);
  std::vector<double> thetas = {0.0, 1e-9, std::numbers::pi - 1e-9};
  for (int i = 0; i < 20; ++i) {
    thetas.push_back(dist(rng));
  }
  double max_err = 0.0;
  for (const double theta : thetas) {
    fill_harmonics(theta, count, s.data(), c.data());
    for (std::size_t j = 1; j <= count; ++j) {
      const double jt = static_cast<double>(j) * theta;
      max_err = std::max(max_err, std::abs(s[j - 1] - std::sin(jt)));
      max_err = std::max(max_err, std::abs(c[j - 1] - std::cos(jt)));
    }
  }
  BOOST_TEST_MESSAGE("max harmonic error: " << max_err);
  BOOST_TEST(max_err < 1e-12);
}

BOOST_AUTO_TEST_CASE(mode_norms_and_orthogonality)
{
  // Non-square domain so that a mix-up of kx and ky would be detected
  const SpectralBasis basis(2.0, 1.5, 3, 3);
  const std::size_t n_modes = basis.n_modes();

  // Midpoint rule; exact (up to rounding) for these trigonometric products
  // because their frequencies are far below the grid's aliasing limit
  constexpr std::size_t n_grid = 64;
  const double hx = basis.lx() / static_cast<double>(n_grid);
  const double hy = basis.ly() / static_cast<double>(n_grid);

  // Gradients of each unit velocity mode at the quadrature points
  std::vector<std::vector<Mat2>> grads(n_modes);
  SpectralBasis::Workspace ws;
  for (std::size_t k = 0; k < n_modes; ++k) {
    std::vector<double> c(n_modes, 0.0);
    c[k] = 1.0;
    grads[k].reserve(n_grid * n_grid);
    for (std::size_t i = 0; i < n_grid; ++i) {
      for (std::size_t j = 0; j < n_grid; ++j) {
        const Vec2 y{
          (static_cast<double>(i) + 0.5) * hx,
          (static_cast<double>(j) + 0.5) * hy};
        grads[k].push_back(basis.evaluate(c, y, ws).dv);
      }
    }
  }

  for (std::size_t m = 1; m <= basis.mx(); ++m) {
    for (std::size_t n = 1; n <= basis.my(); ++n) {
      const std::size_t a = basis.index(m, n);
      const double g = basis.h1_norm_sq(m, n);
      for (std::size_t b = 0; b < n_modes; ++b) {
        double integral = 0.0;
        for (std::size_t p = 0; p < grads[a].size(); ++p) {
          integral += frobenius_inner(grads[a][p], grads[b][p]);
        }
        integral *= hx * hy;
        const double expected = (a == b) ? g : 0.0;
        BOOST_TEST(
          std::abs(integral - expected) < 1e-11 * g,
          "modes " << a << ", " << b << ": integral " << integral
                   << ", expected " << expected);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(derivatives_match_finite_differences)
{
  const SpectralBasis basis(3.0, 2.0, 5, 4);
  std::mt19937 rng(2024);
  const std::vector<double> c = random_coefficients(basis.n_modes(), rng);
  std::uniform_real_distribution<double> ux(0.1, basis.lx() - 0.1);
  std::uniform_real_distribution<double> uy(0.1, basis.ly() - 0.1);
  constexpr double h = 1e-5;

  // Central differences have truncation error O(h^2 k^3) per derivative
  const double tol_v = 1e-6 * coefficient_scale(basis, c, 1);
  const double tol_dv = 1e-6 * coefficient_scale(basis, c, 2);

  SpectralBasis::Workspace ws;
  for (int trial = 0; trial < 50; ++trial) {
    const Vec2 y{ux(rng), uy(rng)};
    const auto f = basis.evaluate(c, y, ws);
    const auto fxp = basis.evaluate(c, {y.x + h, y.y}, ws);
    const auto fxm = basis.evaluate(c, {y.x - h, y.y}, ws);
    const auto fyp = basis.evaluate(c, {y.x, y.y + h}, ws);
    const auto fym = basis.evaluate(c, {y.x, y.y - h}, ws);

    // v = (ds/dy2, -ds/dy1)
    BOOST_TEST(std::abs(f.v.x - (fyp.s - fym.s) / (2 * h)) < tol_v);
    BOOST_TEST(std::abs(f.v.y + (fxp.s - fxm.s) / (2 * h)) < tol_v);

    // dv_ij = dv_i/dy_j, each entry checked independently. In particular,
    // dv22 is compared with a finite difference of v2, not with -dv11.
    BOOST_TEST(std::abs(f.dv.a11 - (fxp.v.x - fxm.v.x) / (2 * h)) < tol_dv);
    BOOST_TEST(std::abs(f.dv.a12 - (fyp.v.x - fym.v.x) / (2 * h)) < tol_dv);
    BOOST_TEST(std::abs(f.dv.a21 - (fxp.v.y - fxm.v.y) / (2 * h)) < tol_dv);
    BOOST_TEST(std::abs(f.dv.a22 - (fyp.v.y - fym.v.y) / (2 * h)) < tol_dv);
  }
}

BOOST_AUTO_TEST_CASE(velocity_is_tangent_to_boundary_and_zero_at_corners)
{
  const SpectralBasis basis(3.0, 2.0, 6, 5);
  std::mt19937 rng(7);
  const std::vector<double> c = random_coefficients(basis.n_modes(), rng);

  // sin(m pi) is not exactly zero in floating point; allow rounding
  const double tol_s = 1e-13 * coefficient_scale(basis, c, 0);
  const double tol_v = 1e-13 * coefficient_scale(basis, c, 1);

  SpectralBasis::Workspace ws;
  constexpr int n_samples = 40;
  for (int i = 0; i <= n_samples; ++i) {
    const double t = static_cast<double>(i) / n_samples;

    // Vertical sides: s = 0 and normal component v1 = 0
    for (const double x : {0.0, basis.lx()}) {
      const auto f = basis.evaluate(c, {x, t * basis.ly()}, ws);
      BOOST_TEST(std::abs(f.s) < tol_s);
      BOOST_TEST(std::abs(f.v.x) < tol_v);
    }

    // Horizontal sides: s = 0 and normal component v2 = 0
    for (const double yy : {0.0, basis.ly()}) {
      const auto f = basis.evaluate(c, {t * basis.lx(), yy}, ws);
      BOOST_TEST(std::abs(f.s) < tol_s);
      BOOST_TEST(std::abs(f.v.y) < tol_v);
    }
  }

  // Corners do not move
  for (const double x : {0.0, basis.lx()}) {
    for (const double yy : {0.0, basis.ly()}) {
      const auto f = basis.evaluate(c, {x, yy}, ws);
      BOOST_TEST(std::abs(f.v.x) < tol_v);
      BOOST_TEST(std::abs(f.v.y) < tol_v);
    }
  }
}
