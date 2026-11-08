// Tests of the splitting integrators of the conformality flow: symplectic
// Euler and Stormer-Verlet (see conformality/flow.hpp)
#include "conformality/flow.hpp"
#include "conformality/spectral_basis.hpp"
#include "conformality_test_maps.hpp"
#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <random>
#include <stdexcept>
#include <vector>

using conformality::advance_points;
using conformality::det;
using conformality::EulerVariant;
using conformality::FlowSettings;
using conformality::FlowStats;
using conformality::frobenius_sq;
using conformality::identity_mat2;
using conformality::Integrator;
using conformality::integrator_from_string;
using conformality::Mat2;
using conformality::SpectralBasis;
using conformality::symplectic_euler_substep;
using conformality::Vec2;
using conformality::velocity_gradient_bound;
using conformality_test::random_coefficients;

namespace
{

FlowSettings with(const Integrator integrator)
{
  FlowSettings settings;
  settings.integrator = integrator;
  return settings;
}

double distance(const Vec2 a, const Vec2 b)
{
  return std::hypot(a.x - b.x, a.y - b.y);
}

double max_point_distance(
  const std::vector<Vec2> &a,
  const std::vector<Vec2> &b)
{
  double m = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    m = std::max(m, distance(a[i], b[i]));
  }
  return m;
}

std::vector<Vec2> random_interior_points(
  const SpectralBasis &basis,
  const int n,
  std::mt19937 &rng)
{
  std::uniform_real_distribution<double> ux(
    0.1 * basis.lx(),
    0.9 * basis.lx());
  std::uniform_real_distribution<double> uy(
    0.1 * basis.ly(),
    0.9 * basis.ly());
  std::vector<Vec2> pts;
  for (int i = 0; i < n; ++i) {
    pts.push_back({ux(rng), uy(rng)});
  }
  return pts;
}

std::vector<Vec2> flowed(
  std::vector<Vec2> pts,
  const SpectralBasis &basis,
  const std::vector<double> &c,
  const double h,
  const std::size_t n,
  const FlowSettings &settings)
{
  const FlowStats stats =
    advance_points(pts, nullptr, basis, c, h, n, settings);
  BOOST_REQUIRE(stats.converged);
  return pts;
}

constexpr Integrator splitting[] = {
  Integrator::symplectic_euler,
  Integrator::stormer_verlet};

}  // namespace

BOOST_AUTO_TEST_CASE(integrator_names_round_trip)
{
  for (const Integrator i :
       {Integrator::implicit_midpoint,
        Integrator::symplectic_euler,
        Integrator::stormer_verlet}) {
    BOOST_TEST((integrator_from_string(to_string(i)) == i));
  }
  BOOST_CHECK_THROW(
    static_cast<void>(integrator_from_string("rk4")),
    std::invalid_argument);
}

// Each side maps to itself, the corners stay fixed, and every tangent has
// unit determinant, also for substeps far larger than the default
BOOST_AUTO_TEST_CASE(splitting_preserves_sides_corners_and_area)
{
  const SpectralBasis basis(2.0, 1.0, 4, 3);
  std::mt19937 rng(31);
  const auto c = random_coefficients(basis.n_modes(), 0.02, rng);
  const double bound = velocity_gradient_bound(basis, c);

  // Points on the sides (and corners), with the coordinates they must keep
  std::vector<Vec2> pts;
  std::vector<bool> keep_x, keep_y;
  auto add = [&](const Vec2 p, const bool kx, const bool ky) {
    pts.push_back(p);
    keep_x.push_back(kx);
    keep_y.push_back(ky);
  };
  for (int i = 0; i <= 20; ++i) {
    const double t = i / 20.0;
    const bool end = i == 0 || i == 20;
    add({0.0, t * basis.ly()}, true, end);
    add({basis.lx(), t * basis.ly()}, true, end);
    add({t * basis.lx(), 0.0}, end, true);
    add({t * basis.lx(), basis.ly()}, end, true);
  }
  for (const Vec2 p : random_interior_points(basis, 20, rng)) {
    add(p, false, false);
  }

  for (const Integrator integrator : splitting) {
    for (const double factor : {0.1, 1.0}) {
      // One substep of size factor / bound, applied 5 times
      std::vector<Vec2> y = pts;
      std::vector<Mat2> tangents(pts.size(), identity_mat2());
      const double h = 5.0 * factor / bound;
      BOOST_REQUIRE(
        advance_points(y, &tangents, basis, c, h, 5, with(integrator))
          .converged);
      const double tol = 1e-13 * basis.lx();
      for (std::size_t i = 0; i < pts.size(); ++i) {
        if (keep_x[i]) {
          BOOST_TEST(std::abs(y[i].x - pts[i].x) < tol);
        }
        if (keep_y[i]) {
          BOOST_TEST(std::abs(y[i].y - pts[i].y) < tol);
        }
        BOOST_TEST(std::abs(det(tangents[i]) - 1.0) < 1e-13);
      }
      // The motion is not trivial
      BOOST_TEST(max_point_distance(y, pts) > 1e-3);
    }
  }
}

// The adjoint variant with -dt inverts a symplectic Euler substep, positions
// and tangents alike
BOOST_AUTO_TEST_CASE(adjoint_euler_inverts_euler)
{
  const SpectralBasis basis(2.0, 1.0, 4, 3);
  std::mt19937 rng(37);
  const auto c = random_coefficients(basis.n_modes(), 0.02, rng);
  const double dt = 0.5 / velocity_gradient_bound(basis, c);
  SpectralBasis::Workspace ws;
  const FlowSettings settings = with(Integrator::symplectic_euler);
  for (const Vec2 p : random_interior_points(basis, 20, rng)) {
    const auto there = symplectic_euler_substep(
      basis,
      c,
      p,
      dt,
      settings,
      ws,
      EulerVariant::implicit_in_y2);
    const auto back = symplectic_euler_substep(
      basis,
      c,
      there.y,
      -dt,
      settings,
      ws,
      EulerVariant::implicit_in_y1);
    BOOST_TEST_REQUIRE(there.converged);
    BOOST_TEST_REQUIRE(back.converged);
    BOOST_TEST(distance(there.y, p) > 1e-4);
    BOOST_TEST(distance(back.y, p) < 1e-13);
    const Mat2 product = back.tangent * there.tangent;
    const Mat2 diff{
      product.a11 - 1.0,
      product.a12,
      product.a21,
      product.a22 - 1.0};
    BOOST_TEST(std::sqrt(frobenius_sq(diff)) < 1e-13);
  }
}

// Global errors over a fixed time: first order for symplectic Euler, second
// order for Stormer-Verlet. The reference is implicit midpoint with many
// substeps.
BOOST_AUTO_TEST_CASE(splitting_converges_at_expected_order)
{
  const SpectralBasis basis(2.0, 1.0, 4, 3);
  std::mt19937 rng(41);
  const auto c = random_coefficients(basis.n_modes(), 0.02, rng);
  const auto start = random_interior_points(basis, 20, rng);
  const double h = 1.0;
  const auto exact =
    flowed(start, basis, c, h, 1024, with(Integrator::implicit_midpoint));
  for (const Integrator integrator : splitting) {
    const double e8 = max_point_distance(
      flowed(start, basis, c, h, 8, with(integrator)),
      exact);
    const double e16 = max_point_distance(
      flowed(start, basis, c, h, 16, with(integrator)),
      exact);
    const double e32 = max_point_distance(
      flowed(start, basis, c, h, 32, with(integrator)),
      exact);
    BOOST_TEST_MESSAGE(
      to_string(integrator)
      << " errors: " << e8 << ", " << e16 << ", " << e32);
    const double expected =
      integrator == Integrator::symplectic_euler ? 2.0 : 4.0;
    BOOST_TEST(e8 / e16 > 0.85 * expected);
    BOOST_TEST(e8 / e16 < 1.15 * expected);
    BOOST_TEST(e16 / e32 > 0.85 * expected);
    BOOST_TEST(e16 / e32 < 1.15 * expected);
  }
}

// The propagated tangent matches central differences of positions; the
// Richardson extrapolation removes their O(delta^2) error
BOOST_AUTO_TEST_CASE(splitting_tangent_matches_finite_differences)
{
  const SpectralBasis basis(2.0, 1.0, 4, 3);
  std::mt19937 rng(43);
  const auto c = random_coefficients(basis.n_modes(), 0.02, rng);
  const auto start = random_interior_points(basis, 10, rng);
  const double h = 1.0;
  const std::size_t n = 4;
  for (const Integrator integrator : splitting) {
    const FlowSettings settings = with(integrator);
    std::vector<Vec2> pts = start;
    std::vector<Mat2> tangents(start.size(), identity_mat2());
    BOOST_REQUIRE(
      advance_points(pts, &tangents, basis, c, h, n, settings).converged);
    auto fd_jacobian = [&](const Vec2 p, const double delta) {
      const auto f = flowed(
        {{p.x + delta, p.y},
         {p.x - delta, p.y},
         {p.x, p.y + delta},
         {p.x, p.y - delta}},
        basis,
        c,
        h,
        n,
        settings);
      return Mat2{
        (f[0].x - f[1].x) / (2 * delta),
        (f[2].x - f[3].x) / (2 * delta),
        (f[0].y - f[1].y) / (2 * delta),
        (f[2].y - f[3].y) / (2 * delta)};
    };
    // Largest error of the Richardson extrapolation from delta and delta / 2
    auto worst = [&](const double delta) {
      double w = 0.0;
      for (std::size_t i = 0; i < start.size(); ++i) {
        const Mat2 j1 = fd_jacobian(start[i], delta);
        const Mat2 j2 = fd_jacobian(start[i], delta / 2);
        const Mat2 d{
          (4 * j2.a11 - j1.a11) / 3 - tangents[i].a11,
          (4 * j2.a12 - j1.a12) / 3 - tangents[i].a12,
          (4 * j2.a21 - j1.a21) / 3 - tangents[i].a21,
          (4 * j2.a22 - j1.a22) / 3 - tangents[i].a22};
        w = std::max(w, std::sqrt(frobenius_sq(d)));
      }
      return w;
    };
    // With an exact tangent, what remains is the O(delta^4) truncation
    // error, which falls by 16 when delta is halved; an error in the tangent
    // would not fall at all. These substeps are large (dt B = 2.5), so the
    // truncation error of symplectic Euler is about 2e-7 at delta = 1e-3.
    const double w1 = worst(1e-3);
    const double w2 = worst(5e-4);
    BOOST_TEST_MESSAGE(
      to_string(integrator) << ": Richardson FD Jacobian error " << w1
                            << " (delta 1e-3), " << w2 << " (delta 5e-4)");
    BOOST_TEST(w2 < 5e-8);
    BOOST_TEST(w1 / w2 > 10.0);
  }
}

// A substep that would fold the rectangle (1 + dt d2s/dy1dy2 <= 0) is
// reported as failed
BOOST_AUTO_TEST_CASE(folding_substep_is_reported)
{
  const SpectralBasis basis(1.0, 1.0, 1, 1);
  const std::vector<double> c = {0.05};
  // d2s/dy1dy2 = c pi^2 cos(pi y1) cos(pi y2), close to c pi^2 near (0, 0)
  const double dt = -3.0 / (0.05 * std::numbers::pi * std::numbers::pi);
  std::vector<Vec2> pts = {{0.05, 0.05}};
  for (const Integrator integrator : splitting) {
    std::vector<Vec2> y = pts;
    BOOST_TEST(!advance_points(y, nullptr, basis, c, dt, 1, with(integrator))
                  .converged);
  }
}
