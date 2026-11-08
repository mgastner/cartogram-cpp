#define BOOST_TEST_MODULE test_conformality_flow
#include "conformality/flow.hpp"
#include "conformality/node_set.hpp"
#include "conformality/objective.hpp"
#include "conformality/spectral_basis.hpp"
#include "conformality_test_maps.hpp"
#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

using conformality::advance_nodes;
using conformality::advance_points;
using conformality::cayley;
using conformality::det;
using conformality::FlowStats;
using conformality::frobenius_sq;
using conformality::h1_descent;
using conformality::identity_mat2;
using conformality::FlowSettings;
using conformality::Mat2;
using conformality::modal_derivatives;
using conformality::NodeSet;
using conformality::objective;
using conformality::SpectralBasis;
using conformality::Vec2;
using conformality::velocity_gradient_bound;
using conformality_test::asymmetric_map_nodes;
using conformality_test::max_abs;
using conformality_test::random_coefficients;

namespace
{

double distance(const Vec2 a, const Vec2 b)
{
  return std::hypot(a.x - b.x, a.y - b.y);
}

// Number of substeps that keeps dt ||Dv|| <= factor everywhere
std::size_t substeps_for(
  const SpectralBasis &basis,
  const std::vector<double> &c,
  double h,
  double factor)
{
  const double steps =
    std::ceil(std::abs(h) * velocity_gradient_bound(basis, c) / factor);
  return std::max<std::size_t>(1, static_cast<std::size_t>(steps));
}

std::vector<Vec2> random_interior_points(
  const SpectralBasis &basis,
  int n,
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

std::vector<Vec2> flowed(
  std::vector<Vec2> pts,
  const SpectralBasis &basis,
  const std::vector<double> &c,
  double h,
  std::size_t n)
{
  const FlowStats stats = advance_points(pts, nullptr, basis, c, h, n);
  BOOST_REQUIRE(stats.converged);
  return pts;
}

}  // namespace

BOOST_AUTO_TEST_CASE(cayley_has_unit_determinant_for_traceless_gradients)
{
  std::mt19937 rng(99);
  std::uniform_real_distribution<double> dist(-3.0, 3.0);
  for (int trial = 0; trial < 100; ++trial) {
    const double a = dist(rng), b = dist(rng), c = dist(rng);
    const Mat2 d{a, b, c, -a};
    const double dt = 0.1 * std::abs(dist(rng));
    BOOST_TEST(std::abs(det(cayley(d, dt)) - 1.0) < 1e-14);
  }

  // Sanity check that the determinant is not 1 by construction
  const Mat2 not_traceless{1.0, 0.0, 0.0, 1.0};
  BOOST_TEST(std::abs(det(cayley(not_traceless, 0.2)) - 1.0) > 0.1);
}

BOOST_AUTO_TEST_CASE(velocity_gradient_bound_holds)
{
  const SpectralBasis basis(3.0, 2.0, 6, 5);
  std::mt19937 rng(5);
  const auto c = random_coefficients(basis.n_modes(), 1.0, rng);
  const double bound = velocity_gradient_bound(basis, c);
  SpectralBasis::Workspace ws;
  for (const Vec2 y : random_interior_points(basis, 200, rng)) {
    const Mat2 dv = basis.evaluate(c, y, ws).dv;
    BOOST_TEST(std::sqrt(frobenius_sq(dv)) <= bound);
  }
}

// Handover test 3: boundary invariance under a single-mode flow
BOOST_AUTO_TEST_CASE(single_mode_flow_preserves_sides_and_corners)
{
  const SpectralBasis basis(2.0, 1.0, 3, 3);
  std::vector<double> c(basis.n_modes(), 0.0);
  c[basis.index(2, 1)] = 0.03;
  const double h = 1.0;
  const std::size_t n = substeps_for(basis, c, h, 0.1);

  std::vector<Vec2> left, right, bottom, top;
  for (int i = 1; i < 20; ++i) {
    const double t = i / 20.0;
    left.push_back({0.0, t * basis.ly()});
    right.push_back({basis.lx(), t * basis.ly()});
    bottom.push_back({t * basis.lx(), 0.0});
    top.push_back({t * basis.lx(), basis.ly()});
  }
  const std::vector<Vec2> corners = {
    {0.0, 0.0},
    {basis.lx(), 0.0},
    {0.0, basis.ly()},
    {basis.lx(), basis.ly()}};

  // The normal velocity vanishes on each side, up to rounding of sin(m pi)
  const double tol = 1e-13 * basis.lx();
  for (const Vec2 p : flowed(left, basis, c, h, n)) {
    BOOST_TEST(std::abs(p.x) < tol);
  }
  for (const Vec2 p : flowed(right, basis, c, h, n)) {
    BOOST_TEST(std::abs(p.x - basis.lx()) < tol);
  }
  for (const Vec2 p : flowed(bottom, basis, c, h, n)) {
    BOOST_TEST(std::abs(p.y) < tol);
  }
  for (const Vec2 p : flowed(top, basis, c, h, n)) {
    BOOST_TEST(std::abs(p.y - basis.ly()) < tol);
  }

  // Points on a side move along it, so the test above is not vacuous
  const auto moved = flowed(left, basis, c, h, n);
  BOOST_TEST(max_point_distance(moved, left) > 1e-3);

  BOOST_TEST(
    max_point_distance(flowed(corners, basis, c, h, n), corners) < tol);
}

// Implicit midpoint has global error O(dt^2): halving dt divides the error
// in positions, and the drift of the Hamiltonian s, by about 4
BOOST_AUTO_TEST_CASE(flow_converges_at_second_order)
{
  const SpectralBasis basis(2.0, 1.0, 4, 3);
  std::mt19937 rng(17);
  const auto c = random_coefficients(basis.n_modes(), 0.02, rng);
  const auto start = random_interior_points(basis, 30, rng);
  const double h = 1.0;

  const auto y8 = flowed(start, basis, c, h, 8);
  const auto y16 = flowed(start, basis, c, h, 16);
  const auto y32 = flowed(start, basis, c, h, 32);
  const auto y64 = flowed(start, basis, c, h, 64);
  const double e8 = max_point_distance(y8, y16);
  const double e16 = max_point_distance(y16, y32);
  const double e32 = max_point_distance(y32, y64);
  BOOST_TEST_MESSAGE(
    "position differences: " << e8 << ", " << e16 << ", " << e32);
  BOOST_TEST(e8 / e16 > 3.5);
  BOOST_TEST(e8 / e16 < 4.5);
  BOOST_TEST(e16 / e32 > 3.5);
  BOOST_TEST(e16 / e32 < 4.5);

  // The stream function is conserved by the exact autonomous flow
  SpectralBasis::Workspace ws;
  auto drift = [&](const std::vector<Vec2> &end) {
    double m = 0.0;
    for (std::size_t i = 0; i < start.size(); ++i) {
      m = std::max(
        m,
        std::abs(
          basis.evaluate(c, end[i], ws).s -
          basis.evaluate(c, start[i], ws).s));
    }
    return m;
  };
  const double d8 = drift(y8), d16 = drift(y16), d32 = drift(y32);
  BOOST_TEST_MESSAGE(
    "stream-function drift: " << d8 << ", " << d16 << ", " << d32);
  BOOST_TEST(d8 / d16 > 3.0);
  BOOST_TEST(d16 / d32 > 3.0);
}

// Compares the propagated tangent with central differences of positions
// advanced through the same substeps. The difference must shrink like
// delta^2, and the Richardson extrapolation (4 J(delta/2) - J(delta)) / 3,
// which cancels the delta^2 term, must match the tangent closely. A mistake
// in Dv or in the tangent update would leave an O(1) discrepancy even if
// det(tangent) = 1.
BOOST_AUTO_TEST_CASE(tangent_matches_finite_differences_of_positions)
{
  const SpectralBasis basis(2.0, 1.0, 4, 3);
  std::mt19937 rng(23);
  const auto c = random_coefficients(basis.n_modes(), 0.02, rng);
  const auto start = random_interior_points(basis, 20, rng);
  const double h = 1.0;
  const std::size_t n = 16;

  std::vector<Vec2> pts = start;
  std::vector<Mat2> tangents(start.size(), identity_mat2());
  BOOST_REQUIRE(advance_points(pts, &tangents, basis, c, h, n).converged);

  auto fd_jacobian = [&](const Vec2 p, double delta) {
    const auto f = flowed(
      {{p.x + delta, p.y},
       {p.x - delta, p.y},
       {p.x, p.y + delta},
       {p.x, p.y - delta}},
      basis,
      c,
      h,
      n);
    return Mat2{
      (f[0].x - f[1].x) / (2 * delta),
      (f[2].x - f[3].x) / (2 * delta),
      (f[0].y - f[1].y) / (2 * delta),
      (f[2].y - f[3].y) / (2 * delta)};
  };
  auto norm_diff = [](const Mat2 &a, const Mat2 &b) {
    const Mat2 d{a.a11 - b.a11, a.a12 - b.a12, a.a21 - b.a21, a.a22 - b.a22};
    return std::sqrt(frobenius_sq(d));
  };

  constexpr double delta = 1e-3;
  double e1 = 0.0, e2 = 0.0, e_richardson = 0.0;
  for (std::size_t i = 0; i < start.size(); ++i) {
    const Mat2 j1 = fd_jacobian(start[i], delta);
    const Mat2 j2 = fd_jacobian(start[i], delta / 2);
    const Mat2 jr{
      (4 * j2.a11 - j1.a11) / 3,
      (4 * j2.a12 - j1.a12) / 3,
      (4 * j2.a21 - j1.a21) / 3,
      (4 * j2.a22 - j1.a22) / 3};
    e1 = std::max(e1, norm_diff(j1, tangents[i]));
    e2 = std::max(e2, norm_diff(j2, tangents[i]));
    e_richardson = std::max(e_richardson, norm_diff(jr, tangents[i]));

    // Area preservation measured from positions alone, to O(delta^2)
    BOOST_TEST(std::abs(det(jr) - 1.0) < 1e-7);
  }
  BOOST_TEST_MESSAGE(
    "FD Jacobian errors: " << e1 << ", " << e2 << "; Richardson "
                           << e_richardson);
  BOOST_TEST(e1 / e2 > 3.5);
  BOOST_TEST(e1 / e2 < 4.5);
  BOOST_TEST(e_richardson < 1e-7);

  // Tautological but cheap: the propagated tangent has unit determinant
  for (const Mat2 &t : tangents) {
    BOOST_TEST(std::abs(det(t) - 1.0) < 1e-13);
  }

  // The motion is not trivial
  BOOST_TEST(max_point_distance(pts, start) > 1e-2);
}

BOOST_AUTO_TEST_CASE(node_jacobians_are_left_multiplied_by_the_tangent)
{
  const double lx = 2.0, ly = 1.0;
  const NodeSet before = asymmetric_map_nodes(lx, ly, 16, 8);
  const SpectralBasis basis(lx, ly, 4, 3);
  std::mt19937 rng(31);
  const auto c = random_coefficients(basis.n_modes(), 0.02, rng);

  NodeSet after = before;
  BOOST_REQUIRE(advance_nodes(after, basis, c, 0.5, 8).converged);
  std::vector<Vec2> pts = before.y;
  std::vector<Mat2> tangents(pts.size(), identity_mat2());
  BOOST_REQUIRE(advance_points(pts, &tangents, basis, c, 0.5, 8).converged);
  for (std::size_t j = 0; j < before.size(); ++j) {
    BOOST_TEST(distance(after.y[j], pts[j]) < 1e-15);
    const Mat2 expected = tangents[j] * before.jac[j];
    const Mat2 diff{
      after.jac[j].a11 - expected.a11,
      after.jac[j].a12 - expected.a12,
      after.jac[j].a21 - expected.a21,
      after.jac[j].a22 - expected.a22};
    BOOST_TEST(std::sqrt(frobenius_sq(diff)) < 1e-13);
    BOOST_TEST(std::abs(det(after.jac[j]) / before.g0[j] - 1.0) < 1e-13);
  }
}

BOOST_AUTO_TEST_CASE(non_convergent_solve_is_reported)
{
  const SpectralBasis basis(1.0, 1.0, 3, 3);
  std::mt19937 rng(3);
  const auto c = random_coefficients(basis.n_modes(), 1.0, rng);

  // One huge substep: the fixed-point map is far from a contraction
  FlowSettings settings;
  settings.max_iterations = 5;
  std::vector<Vec2> pts = {{0.3, 0.4}, {0.6, 0.7}};
  const FlowStats stats =
    advance_points(pts, nullptr, basis, c, 50.0, 1, settings);
  BOOST_TEST(!stats.converged);
}

// Handover test 4: F_mn against central differences of the objective along
// the actual flow, E0(A_eta o T) - E0(A_-eta o T) over 2 eta
BOOST_AUTO_TEST_CASE(modal_derivatives_match_flow_finite_differences)
{
  const double lx = 2.0, ly = 1.0;
  const NodeSet nodes = asymmetric_map_nodes(lx, ly, 64, 32);
  const SpectralBasis basis(lx, ly, 5, 4);
  const std::vector<double> f = modal_derivatives(nodes, basis);
  const double scale = max_abs(f);
  BOOST_TEST(scale > 1e-3);

  auto flowed_energy = [&](const std::vector<double> &c, double h) {
    NodeSet copy = nodes;
    BOOST_REQUIRE(advance_nodes(copy, basis, c, h, 2).converged);
    return objective(copy);
  };

  for (std::size_t m = 1; m <= basis.mx(); ++m) {
    for (std::size_t n = 1; n <= basis.my(); ++n) {
      std::vector<double> c(basis.n_modes(), 0.0);
      const std::size_t k = basis.index(m, n);
      c[k] = 1.0;
      const double k_sq =
        basis.kx(m) * basis.kx(m) + basis.ky(n) * basis.ky(n);
      const double eta = 1e-4 / k_sq;
      const double fd =
        (flowed_energy(c, eta) - flowed_energy(c, -eta)) / (2 * eta);
      BOOST_TEST(
        std::abs(fd - f[k]) < 1e-6 * scale,
        "mode (" << m << ", " << n << "): FD " << fd << ", F " << f[k]);
    }
  }

  // The descent velocity has slope -R along its own flow
  const auto descent = h1_descent(basis, f);
  const double eta = 1e-4 / velocity_gradient_bound(basis, descent.coeffs);
  const double fd = (flowed_energy(descent.coeffs, eta) -
                     flowed_energy(descent.coeffs, -eta)) /
                    (2 * eta);
  BOOST_TEST_MESSAGE("R = " << descent.r << ", flow FD slope = " << fd);
  BOOST_TEST(std::abs(fd + descent.r) < 1e-6 * descent.r);
}
