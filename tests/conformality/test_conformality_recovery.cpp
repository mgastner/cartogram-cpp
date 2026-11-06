#define BOOST_TEST_MODULE test_conformality_recovery
#include "conformality/flow.hpp"
#include "conformality/node_set.hpp"
#include "conformality/objective.hpp"
#include "conformality/optimizer.hpp"
#include "conformality/spectral_basis.hpp"
#include "conformality_test_maps.hpp"
#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <cstddef>
#include <vector>

using conformality::advance_nodes;
using conformality::det;
using conformality::make_cell_centre_nodes;
using conformality::max_det_drift;
using conformality::NodeSet;
using conformality::objective;
using conformality::OptimizationResult;
using conformality::optimize;
using conformality::OptimizerSettings;
using conformality::Outcome;
using conformality::replay;
using conformality::SpectralBasis;
using conformality::to_string;
using conformality::Vec2;
using conformality_test::asymmetric_map_nodes;
using conformality_test::shear_map_nodes;

namespace
{

double max_distance(const std::vector<Vec2> &a, const std::vector<Vec2> &b)
{
  double m = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    m = std::max(m, std::hypot(a[i].x - b[i].x, a[i].y - b[i].y));
  }
  return m;
}

// Stream coefficients of a small three-mode area-preserving flow B
std::vector<double> perturbation_coefficients(const SpectralBasis &basis)
{
  std::vector<double> c(basis.n_modes(), 0.0);
  c[basis.index(1, 1)] = 0.02;
  c[basis.index(2, 1)] = -0.012;
  c[basis.index(1, 2)] = 0.008;
  return c;
}

void log_summary(const OptimizationResult &result)
{
  BOOST_TEST_MESSAGE(
    "outcome " << to_string(result.outcome) << " (" << result.reason << "), "
               << result.steps.size() << " accepted steps, E "
               << result.initial_energy << " -> " << result.final_energy
               << ", sqrt(R/R0) "
               << std::sqrt(result.final_r / result.initial_r));
  const std::size_t n = result.history.size();
  for (std::size_t i = 0; i < n; ++i) {
    if (i < 5 || i % 25 == 0 || i + 1 == n) {
      const auto &rec = result.history[i];
      BOOST_TEST_MESSAGE(
        "  it " << rec.iteration << ": E " << rec.energy << ", R " << rec.r
                << ", h " << rec.h << ", substeps " << rec.n_substeps
                << ", backtracks " << rec.backtracks << ", E/E_prev "
                << rec.energy_after / rec.energy);
    }
  }
}

// Accepted steps satisfy Armijo, energies decrease monotonically, and the
// frozen Jacobian determinant is preserved
void check_history(const OptimizationResult &result, double det_tol)
{
  double previous = result.initial_energy;
  for (const auto &rec : result.history) {
    if (rec.h > 0.0) {
      BOOST_TEST(rec.armijo_margin >= 0.0);
      BOOST_TEST(rec.energy_after < previous);
      BOOST_TEST(rec.max_det_drift <= det_tol);
      previous = rec.energy_after;
    }
  }
}

}  // namespace

// Revised handover test 5: T0 = B is area preserving (g0 = 1), so the
// unique admissible minimizer is the identity with E0 = 0. B uses three
// modes of a 4 x 4 basis; the optimizer uses 8 x 8 modes. The optimizer
// only reaches stationarity in its retained modes, so the remaining
// residual shrinks as modes are added (about 1e-8 at 4 x 4, 3e-13 at
// 8 x 8 when the test was written).
BOOST_AUTO_TEST_CASE(recovers_identity_from_area_preserving_map)
{
  const double lx = 2.0, ly = 1.0;
  const SpectralBasis b_basis(lx, ly, 4, 4);
  const SpectralBasis basis(lx, ly, 8, 8);
  NodeSet nodes = make_cell_centre_nodes(lx, ly, 64, 32);
  BOOST_REQUIRE(
    advance_nodes(nodes, b_basis, perturbation_coefficients(b_basis), 1.0, 64)
      .converged);
  for (std::size_t j = 0; j < nodes.size(); ++j) {
    nodes.g0[j] = det(nodes.jac[j]);
  }
  const NodeSet start = nodes;
  const double e_start = objective(start);
  const double d_start = max_distance(start.y, start.x);
  BOOST_TEST_MESSAGE(
    "start: E0 = " << e_start << ", max |T(x) - x| = " << d_start);

  OptimizerSettings settings;
  settings.gradient_tol = 1e-6;
  const OptimizationResult result = optimize(nodes, basis, settings);
  log_summary(result);
  check_history(result, settings.det_drift_tol);

  const double d_final = max_distance(nodes.y, nodes.x);
  BOOST_TEST_MESSAGE("final: max |T(x) - x| = " << d_final);
  BOOST_TEST(result.outcome == Outcome::converged);
  BOOST_TEST(result.final_energy <= 1e-9 * e_start);
  BOOST_TEST(d_final <= 1e-5 * d_start);
}

namespace
{

struct ShearRecovery {
  OptimizationResult result;
  double e_p{0.0};
  double e_start{0.0};
  double d_start{0.0};
  double d_final{0.0};
};

// T0 = B o P with B from three modes of a 4 x 4 basis, optimized with
// modes x modes retained modes
ShearRecovery recover_shear_map(std::size_t modes)
{
  constexpr double eps = 0.3;
  const SpectralBasis b_basis(1.0, 1.0, 4, 4);
  const SpectralBasis basis(1.0, 1.0, modes, modes);
  const NodeSet p = shear_map_nodes(eps, 64);
  ShearRecovery out;
  out.e_p = objective(p);
  BOOST_TEST(
    std::abs(out.e_p - (1.0 / std::sqrt(1.0 - eps * eps) - 1.0)) < 1e-13);

  NodeSet nodes = p;
  BOOST_REQUIRE(
    advance_nodes(nodes, b_basis, perturbation_coefficients(b_basis), 1.0, 64)
      .converged);
  out.e_start = objective(nodes);
  out.d_start = max_distance(nodes.y, p.y);
  BOOST_REQUIRE(out.e_start > out.e_p);

  OptimizerSettings settings;
  settings.gradient_tol = 1e-6;
  out.result = optimize(nodes, basis, settings);
  log_summary(out.result);
  check_history(out.result, settings.det_drift_tol);
  out.d_final = max_distance(nodes.y, p.y);
  BOOST_TEST_MESSAGE(
    modes << " x " << modes
          << " modes: E0 - E0[P] = " << out.result.final_energy - out.e_p
          << ", max |T(x) - P(x)| = " << out.d_final << " (start "
          << out.d_start << ")");
  BOOST_TEST(out.result.outcome == Outcome::converged);
  return out;
}

}  // namespace

// Revised handover test 6: T0 = B o P with P the shear map, which is a
// strict local minimizer for its own frozen Jacobian. The optimizer must
// return to P, with E0[P] = (1 - eps^2)^(-1/2) - 1, up to a residual that
// is set by the retained modes and shrinks as modes are added.
BOOST_AUTO_TEST_CASE(recovers_shear_map_with_non_uniform_jacobian)
{
  const ShearRecovery coarse = recover_shear_map(4);
  const ShearRecovery fine = recover_shear_map(8);

  // Two-sided: the energy returns to E0[P] and the map returns to P
  const double gap_coarse = coarse.result.final_energy - coarse.e_p;
  const double gap_fine = fine.result.final_energy - fine.e_p;
  BOOST_TEST(std::abs(gap_fine) <= 1e-5 * (fine.e_start - fine.e_p));
  BOOST_TEST(fine.d_final <= 1e-3 * fine.d_start);

  // Refinement in the retained modes reduces the residual
  BOOST_TEST(std::abs(gap_fine) * 10.0 < std::abs(gap_coarse));
  BOOST_TEST(fine.d_final * 10.0 < coarse.d_final);
}

// Replaying the logged steps on a copy of the starting nodes reproduces the
// optimized nodes bit for bit, and the optimizer leaves g0 untouched
BOOST_AUTO_TEST_CASE(replay_reproduces_optimized_nodes)
{
  const double lx = 2.0, ly = 1.0;
  const SpectralBasis basis(lx, ly, 4, 3);
  const NodeSet start = asymmetric_map_nodes(lx, ly, 32, 16);
  NodeSet nodes = start;
  OptimizerSettings settings;
  settings.max_iterations = 15;
  const OptimizationResult result = optimize(nodes, basis, settings);
  log_summary(result);
  check_history(result, settings.det_drift_tol);
  BOOST_TEST(result.outcome != Outcome::numerical_failure);
  BOOST_REQUIRE(!result.steps.empty());
  BOOST_TEST(result.final_energy < result.initial_energy);

  NodeSet replayed = start;
  BOOST_REQUIRE(replay(replayed, basis, result.steps).converged);
  double max_diff = 0.0;
  for (std::size_t j = 0; j < nodes.size(); ++j) {
    max_diff = std::max(
      {max_diff,
       std::abs(replayed.y[j].x - nodes.y[j].x),
       std::abs(replayed.y[j].y - nodes.y[j].y),
       std::abs(replayed.jac[j].a11 - nodes.jac[j].a11),
       std::abs(replayed.jac[j].a12 - nodes.jac[j].a12),
       std::abs(replayed.jac[j].a21 - nodes.jac[j].a21),
       std::abs(replayed.jac[j].a22 - nodes.jac[j].a22)});
  }
  BOOST_TEST(max_diff <= 0.0);
  BOOST_TEST(max_det_drift(nodes) <= settings.det_drift_tol);
  for (std::size_t j = 0; j < nodes.size(); ++j) {
    BOOST_TEST(std::abs(nodes.g0[j] - start.g0[j]) <= 0.0);
  }
}
