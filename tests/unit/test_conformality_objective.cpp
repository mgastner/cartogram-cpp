#define BOOST_TEST_MODULE test_conformality_objective
#include "conformality/node_set.hpp"
#include "conformality/objective.hpp"
#include "conformality/spectral_basis.hpp"
#include "conformality_test_maps.hpp"
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

using conformality::det;
using conformality::excess_mips;
using conformality::h1_descent;
using conformality::identity_mat2;
using conformality::make_cell_centre_nodes;
using conformality::Mat2;
using conformality::modal_derivatives;
using conformality::NodeSet;
using conformality::objective;
using conformality::SpectralBasis;
using conformality_test::asymmetric_map_nodes;
using conformality_test::max_abs;
using conformality_test::shear_map_nodes;

namespace
{

// Nodes whose Jacobians are perturbed by J -> (I + eta Dv(y)) J, where Dv is
// the velocity gradient of the stream function with coefficients c
NodeSet perturbed_nodes(
  const NodeSet &nodes,
  const SpectralBasis &basis,
  const std::vector<double> &c,
  double eta)
{
  NodeSet out = nodes;
  SpectralBasis::Workspace ws;
  for (std::size_t j = 0; j < nodes.size(); ++j) {
    const Mat2 dv = basis.evaluate(c, nodes.y[j], ws).dv;
    const Mat2 step{
      1.0 + eta * dv.a11,
      eta * dv.a12,
      eta * dv.a21,
      1.0 + eta * dv.a22};
    out.jac[j] = step * nodes.jac[j];
  }
  return out;
}

}  // namespace

BOOST_AUTO_TEST_CASE(excess_mips_values)
{
  BOOST_TEST(std::abs(excess_mips(identity_mat2())) < 1e-15);

  // A scaled rotation is conformal
  const double t = 0.7;
  const Mat2 similarity{
    3.0 * std::cos(t),
    -3.0 * std::sin(t),
    3.0 * std::sin(t),
    3.0 * std::cos(t)};
  BOOST_TEST(std::abs(excess_mips(similarity)) < 1e-14);

  // Singular values 2 and 1: K + 1/K - 2 = 0.5, independent of scale
  BOOST_TEST(std::abs(excess_mips({2.0, 0.0, 0.0, 1.0}) - 0.5) < 1e-15);
  BOOST_TEST(std::abs(excess_mips({10.0, 0.0, 0.0, 5.0}) - 0.5) < 1e-15);

  // Nearly conformal J: diag(1 + e, 1) has K + 1/K - 2 = e^2 / (1 + e)
  // exactly. The value must keep full relative accuracy, which the naive
  // form ||J||^2 / det J - 2 loses to cancellation.
  const double e = 1e-7;
  const double near = excess_mips({1.0 + e, 0.0, 0.0, 1.0});
  BOOST_TEST(std::abs(near / (e * e / (1.0 + e)) - 1.0) < 1e-8);

  // A general matrix agrees with the defining formula
  const Mat2 general{1.3, -0.4, 0.7, 0.9};
  BOOST_TEST(
    std::abs(
      excess_mips(general) -
      (conformality::frobenius_sq(general) / det(general) - 2.0)) < 1e-14);
}

BOOST_AUTO_TEST_CASE(distortion_measures_match_singular_values)
{
  using conformality::log_axis_ratio;
  using conformality::max_angular_distortion;

  // Rotation times diag(3, 1.5) times rotation: a / b = 2
  const double t = 0.4, u = -1.1;
  const Mat2 r1{std::cos(t), -std::sin(t), std::sin(t), std::cos(t)};
  const Mat2 r2{std::cos(u), -std::sin(u), std::sin(u), std::cos(u)};
  const Mat2 j = r1 * Mat2{3.0, 0.0, 0.0, 1.5} * r2;
  BOOST_TEST(std::abs(log_axis_ratio(j) - std::log(2.0)) < 1e-14);
  BOOST_TEST(
    std::abs(max_angular_distortion(j) - 2.0 * std::asin(1.0 / 3.0)) < 1e-14);

  // Conformal J has no distortion; nearly conformal J keeps relative
  // accuracy: ln(1 + e) for diag(1 + e, 1)
  BOOST_TEST(std::abs(log_axis_ratio(r1)) < 1e-15);
  const double e = 1e-9;
  BOOST_TEST(
    std::abs(log_axis_ratio({1.0 + e, 0.0, 0.0, 1.0}) / std::log1p(e) - 1.0) <
    1e-6);
}

BOOST_AUTO_TEST_CASE(identity_map_has_zero_energy_and_derivatives)
{
  const NodeSet nodes = make_cell_centre_nodes(2.0, 1.0, 40, 20);
  const SpectralBasis basis(2.0, 1.0, 6, 5);
  BOOST_TEST(std::abs(objective(nodes)) < 1e-15);
  BOOST_TEST(max_abs(modal_derivatives(nodes, basis)) < 1e-15);
}

BOOST_AUTO_TEST_CASE(objective_rejects_inverted_or_degenerate_nodes)
{
  NodeSet nodes = make_cell_centre_nodes(1.0, 1.0, 4, 4);
  nodes.jac[5] = {1.0, 0.0, 0.0, -1.0};
  BOOST_TEST(std::isinf(objective(nodes)));
  const SpectralBasis basis(1.0, 1.0, 2, 2);
  BOOST_CHECK_THROW(
    static_cast<void>(modal_derivatives(nodes, basis)),
    std::domain_error);
  nodes.jac[5] = {1.0, 1.0, 1.0, 1.0};
  BOOST_TEST(std::isinf(objective(nodes)));
}

// Handover test 2. The map is separable, so for alpha = 0 every modal
// derivative vanishes; the energy has the closed form (1 - eps^2)^(-1/2) - 1.
BOOST_AUTO_TEST_CASE(shear_map_matches_closed_form_and_is_stationary)
{
  const SpectralBasis basis(1.0, 1.0, 8, 8);
  for (const double eps : {0.3, 0.7}) {
    const NodeSet nodes = shear_map_nodes(eps, 128);
    const double exact = 1.0 / std::sqrt(1.0 - eps * eps) - 1.0;
    const double e = objective(nodes);
    BOOST_TEST_MESSAGE(
      "eps = " << eps << ": E0 = " << e << ", exact " << exact);
    BOOST_TEST(std::abs(e - exact) < 1e-13);

    // Scale of individual terms in F: kmax^2 times the energy density
    const double kmax = basis.kx(basis.mx());
    const double f_scale = kmax * kmax * (1.0 + exact);
    BOOST_TEST(max_abs(modal_derivatives(nodes, basis)) < 1e-12 * f_scale);
  }
}

// Compares the tabulated assembly of F_mn with a direct sum over nodes that
// uses SpectralBasis::evaluate for each unit mode
BOOST_AUTO_TEST_CASE(modal_derivatives_match_reference_summation)
{
  const double lx = 2.0, ly = 1.0;
  const NodeSet nodes = asymmetric_map_nodes(lx, ly, 96, 48);
  for (const auto &jac : nodes.jac) {
    BOOST_REQUIRE(det(jac) > 0.0);
  }
  const SpectralBasis basis(lx, ly, 6, 5);
  const std::vector<double> f = modal_derivatives(nodes, basis);

  double z = 0.0;
  for (const double w : nodes.w) {
    z += w;
  }
  SpectralBasis::Workspace ws;
  std::vector<double> f_ref(basis.n_modes(), 0.0);
  for (std::size_t k = 0; k < basis.n_modes(); ++k) {
    std::vector<double> c(basis.n_modes(), 0.0);
    c[k] = 1.0;
    for (std::size_t j = 0; j < nodes.size(); ++j) {
      const Mat2 &jac = nodes.jac[j];
      const Mat2 dv = basis.evaluate(c, nodes.y[j], ws).dv;
      const double a11 = jac.a11 * jac.a11 + jac.a12 * jac.a12;
      const double a22 = jac.a21 * jac.a21 + jac.a22 * jac.a22;
      const double a12 = jac.a11 * jac.a21 + jac.a12 * jac.a22;
      const double contraction =
        a11 * dv.a11 + a12 * dv.a12 + a12 * dv.a21 + a22 * dv.a22;
      f_ref[k] += 2.0 * nodes.w[j] * contraction / (z * det(jac));
    }
  }
  const double scale = max_abs(f_ref);
  BOOST_TEST_MESSAGE("max |F_mn| = " << scale);
  BOOST_TEST(scale > 1e-3);
  for (std::size_t k = 0; k < f.size(); ++k) {
    BOOST_TEST(
      std::abs(f[k] - f_ref[k]) < 1e-12 * scale,
      "mode " << k << ": " << f[k] << " vs " << f_ref[k]);
  }
}

// F_mn is the derivative of the discrete objective along J -> (I + eta Dv) J.
// Because Dv is traceless, this has the same first derivative at eta = 0 as
// the flow of v, so the check is independent of any time integrator.
BOOST_AUTO_TEST_CASE(modal_derivatives_match_finite_differences)
{
  const double lx = 2.0, ly = 1.0;
  const NodeSet nodes = asymmetric_map_nodes(lx, ly, 96, 48);
  const SpectralBasis basis(lx, ly, 6, 5);
  const std::vector<double> f = modal_derivatives(nodes, basis);
  const double scale = max_abs(f);

  for (std::size_t m = 1; m <= basis.mx(); ++m) {
    for (std::size_t n = 1; n <= basis.my(); ++n) {
      std::vector<double> c(basis.n_modes(), 0.0);
      const std::size_t k = basis.index(m, n);
      c[k] = 1.0;

      // |Dv_mn| grows like kx^2 + ky^2; scale the step so that eta |Dv| is
      // about 1e-4 for every mode, keeping the O(eta^2) error uniform
      const double k_sq =
        basis.kx(m) * basis.kx(m) + basis.ky(n) * basis.ky(n);
      const double eta = 1e-4 / k_sq;
      const double ep = objective(perturbed_nodes(nodes, basis, c, eta));
      const double em = objective(perturbed_nodes(nodes, basis, c, -eta));
      const double fd = (ep - em) / (2 * eta);
      BOOST_TEST(
        std::abs(fd - f[k]) < 1e-6 * scale,
        "mode (" << m << ", " << n << "): FD " << fd << ", F " << f[k]);
    }
  }
}

// The descent velocity has initial slope -R, both as the sum c . F and as a
// finite difference of the objective
BOOST_AUTO_TEST_CASE(descent_direction_has_slope_minus_r)
{
  const double lx = 2.0, ly = 1.0;
  const NodeSet nodes = asymmetric_map_nodes(lx, ly, 96, 48);
  const SpectralBasis basis(lx, ly, 6, 5);
  const std::vector<double> f = modal_derivatives(nodes, basis);
  const auto descent = h1_descent(basis, f);
  BOOST_TEST(descent.r > 0.0);

  double slope = 0.0;
  for (std::size_t m = 1; m <= basis.mx(); ++m) {
    for (std::size_t n = 1; n <= basis.my(); ++n) {
      const std::size_t k = basis.index(m, n);
      const double g = basis.h1_norm_sq(m, n);
      BOOST_TEST(
        std::abs(descent.coeffs[k] * g + f[k]) <= 1e-14 * std::abs(f[k]));
      slope += descent.coeffs[k] * f[k];
    }
  }
  BOOST_TEST(std::abs(slope + descent.r) < 1e-13 * descent.r);

  constexpr double eta = 1e-4;
  const double ep =
    objective(perturbed_nodes(nodes, basis, descent.coeffs, eta));
  const double em =
    objective(perturbed_nodes(nodes, basis, descent.coeffs, -eta));
  const double fd = (ep - em) / (2 * eta);
  BOOST_TEST_MESSAGE("R = " << descent.r << ", FD slope = " << fd);
  BOOST_TEST(std::abs(fd + descent.r) < 1e-6 * descent.r);
}
