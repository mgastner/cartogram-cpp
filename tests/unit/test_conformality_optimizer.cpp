#define BOOST_TEST_MODULE test_conformality_optimizer
#include "conformality/flow.hpp"
#include "conformality/mesh_check.hpp"
#include "conformality/node_set.hpp"
#include "conformality/objective.hpp"
#include "conformality/optimizer.hpp"
#include "conformality/parallel.hpp"
#include "conformality/spectral_basis.hpp"
#include "conformality_test_maps.hpp"
#include <algorithm>
#include <bit>
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using conformality::advance_nodes;
using conformality::max_threads;
using conformality::modal_derivatives;
using conformality::NodeSet;
using conformality::openmp_enabled;
using conformality::optimize;
using conformality::OptimizerSettings;
using conformality::Outcome;
using conformality::set_max_threads;
using conformality::SpectralBasis;
using conformality::velocity_gradient_bound;
using conformality_test::asymmetric_map_nodes;
using conformality_test::random_coefficients;

namespace
{

constexpr double lx = 2.0;
constexpr double ly = 1.0;

// Whether two sequences of doubles agree bit for bit
bool same_bits(const std::vector<double> &a, const std::vector<double> &b)
{
  if (a.size() != b.size()) {
    return false;
  }
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (
      std::bit_cast<std::uint64_t>(a[i]) !=
      std::bit_cast<std::uint64_t>(b[i])) {
      return false;
    }
  }
  return true;
}

// Restores the thread count when a test ends
class ThreadCountGuard
{
public:
  ThreadCountGuard() : saved_(max_threads()) {}
  ~ThreadCountGuard()
  {
    set_max_threads(saved_);
  }
  ThreadCountGuard(const ThreadCountGuard &) = delete;
  ThreadCountGuard &operator=(const ThreadCountGuard &) = delete;

private:
  int saved_;
};

}  // namespace

// Every trial, and hence every accepted step, satisfies h B <= max_strain.
// With a tiny cap, the first trial is the capped one and is accepted.
BOOST_AUTO_TEST_CASE(strain_cap_bounds_every_step)
{
  NodeSet nodes = asymmetric_map_nodes(lx, ly, 32, 16);
  const SpectralBasis basis(lx, ly, 4, 2);
  OptimizerSettings settings;
  settings.max_strain = 1e-3;
  settings.max_iterations = 4;
  settings.stall_window = 0;
  const auto result = optimize(nodes, basis, settings);
  BOOST_TEST_REQUIRE(result.history.size() == 4u);
  BOOST_TEST(result.history[0].backtracks == 0u);
  BOOST_TEST(
    result.history[0].strain == settings.max_strain,
    boost::test_tools::tolerance(1e-12));
  for (std::size_t k = 0; k < result.history.size(); ++k) {
    const auto &r = result.history[k];
    const double bound =
      velocity_gradient_bound(basis, result.steps[k].coeffs);
    BOOST_TEST(r.strain <= settings.max_strain * (1.0 + 1e-12));
    BOOST_TEST(r.strain == r.h * bound, boost::test_tools::tolerance(1e-14));
    BOOST_TEST(r.energy_after < r.energy);
  }
}

// The stall test stops after the first accepted step k with
// E_(k - w) - E_k <= stall_tol * E_k, where w is the window. A reference
// run without the test gives E_0, ..., E_6; a tolerance between the
// relative window decreases q_k = (E_(k-3) - E_k) / E_k then fixes the step
// at which a run with w = 3 must stop.
BOOST_AUTO_TEST_CASE(stall_window_stops_slow_tail)
{
  const NodeSet start = asymmetric_map_nodes(lx, ly, 32, 16);
  const SpectralBasis basis(lx, ly, 4, 2);

  // A window of 0 disables the test
  OptimizerSettings settings;
  settings.stall_window = 0;
  settings.energy_tol = 0.0;
  settings.max_iterations = 6;
  NodeSet nodes = start;
  const auto reference = optimize(nodes, basis, settings);
  BOOST_TEST_REQUIRE(reference.outcome == Outcome::budget_exhausted);
  BOOST_TEST_REQUIRE(reference.history.size() == 6u);
  std::vector<double> e;
  for (const auto &r : reference.history) {
    e.push_back(r.energy);
  }
  e.push_back(reference.final_energy);
  auto q = [&](const std::size_t k) {
    return (e[k - 3] - e[k]) / e[k];
  };
  BOOST_TEST_MESSAGE(
    "window decreases q3 " << q(3) << ", q4 " << q(4) << ", q5 " << q(5));
  const double earlier = std::min(q(3), q(4));
  BOOST_TEST_REQUIRE(q(5) < earlier);

  settings.stall_window = 3;
  settings.stall_tol = std::sqrt(q(5) * earlier);
  nodes = start;
  const auto stalled = optimize(nodes, basis, settings);
  BOOST_TEST(stalled.outcome == Outcome::converged);
  BOOST_TEST(stalled.history.size() == 5u);
  BOOST_TEST(stalled.reason.find("last 3 iterations") != std::string::npos);
  BOOST_TEST(same_bits({stalled.final_energy}, {e[5]}));
}

// Without steps, every cell keeps its area up to rounding. A single smooth
// area-preserving step bends the cell edges, so straight-edged image cells
// have an O((h / s)^2) area error for cells of size h with s samples per
// edge: in the asymptotic regime it falls by about 4 per halving of h and
// per doubling of s.
BOOST_AUTO_TEST_CASE(mesh_check_converges_at_second_order)
{
  const SpectralBasis basis(lx, ly, 2, 1);
  const auto identity =
    conformality::check_mesh_areas(basis, {}, lx, ly, 4, 4);
  BOOST_TEST_REQUIRE(identity.size() == 5u);
  for (const auto &c : identity) {
    BOOST_TEST(c.max_rel_error < 1e-14);
  }

  // A generic smooth flow: no cell edge lies on a nodal line
  const SpectralBasis generic(lx, ly, 3, 2);
  std::mt19937 rng(11);
  conformality::AcceptedStep step;
  step.coeffs = random_coefficients(generic.n_modes(), 1e-2, rng);
  step.h = 1.0;
  step.n_substeps = 4;
  const auto checks =
    conformality::check_mesh_areas(generic, {step}, lx, ly, 4, 8);
  for (const auto &c : checks) {
    BOOST_TEST_MESSAGE(
      c.cells_long << " cells, " << c.samples_per_edge
                   << " samples per edge: " << c.max_rel_error);
  }
  BOOST_TEST_REQUIRE(checks.size() == 7u);

  // Finer meshes with one sample per edge (4, 8, 16 and 32 cells): the
  // largest error moves between cells until they are small compared with
  // the length scale of the bending, so only the finest pair is asymptotic
  for (std::size_t k = 1; k < 4; ++k) {
    BOOST_TEST(checks[k].cells_long == 2 * checks[k - 1].cells_long);
    BOOST_TEST(checks[k].max_rel_error < checks[k - 1].max_rel_error);
  }
  const double finest = checks[2].max_rel_error / checks[3].max_rel_error;
  BOOST_TEST(finest > 3.0);
  BOOST_TEST(finest < 5.0);

  // Four cells with 1, 2, 4 and 8 samples per edge
  const std::vector<double> samples = {
    checks[0].max_rel_error,
    checks[4].max_rel_error,
    checks[5].max_rel_error,
    checks[6].max_rel_error};
  BOOST_TEST(checks[6].samples_per_edge == 8u);
  for (std::size_t k = 1; k < samples.size(); ++k) {
    const double ratio = samples[k - 1] / samples[k];
    BOOST_TEST(ratio > 3.0);
    BOOST_TEST(ratio < 5.0);
  }

  BOOST_CHECK_THROW(
    static_cast<void>(
      conformality::check_mesh_areas(generic, {step}, lx, ly, 4, 3)),
    std::invalid_argument);
}

// A seed gives the same random step on every call. The step meets its strain
// bound, uses the optimizer's substep rule and preserves area.
BOOST_AUTO_TEST_CASE(random_step_is_reproducible_and_area_preserving)
{
  const SpectralBasis basis(lx, ly, 6, 4);
  const auto a = conformality::random_step(basis, 2.0, 7);
  const auto b = conformality::random_step(basis, 2.0, 7);
  const auto c = conformality::random_step(basis, 2.0, 8);
  const bool same =
    (a.coeffs == b.coeffs) &&
    (std::bit_cast<std::uint64_t>(a.h) == std::bit_cast<std::uint64_t>(b.h));
  const bool different = (a.coeffs != c.coeffs);
  BOOST_TEST(same);
  BOOST_TEST(different);
  BOOST_TEST(
    std::abs(a.h * velocity_gradient_bound(basis, a.coeffs) - 2.0) < 1e-12);
  BOOST_TEST(a.n_substeps == 20u);

  NodeSet nodes = conformality::make_cell_centre_nodes(lx, ly, 16, 8);
  const NodeSet start = nodes;
  const auto stats = conformality::replay(nodes, basis, {a});
  BOOST_TEST(stats.converged);
  BOOST_TEST(conformality::max_det_drift(nodes) < 1e-12);
  double moved = 0.0;
  for (std::size_t j = 0; j < nodes.size(); ++j) {
    moved = std::max(
      moved,
      std::hypot(nodes.y[j].x - start.y[j].x, nodes.y[j].y - start.y[j].y));
  }
  BOOST_TEST(moved > 1e-3);

  BOOST_CHECK_THROW(
    static_cast<void>(conformality::random_step(basis, 0.0, 1)),
    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(invalid_settings_are_rejected)
{
  NodeSet nodes = asymmetric_map_nodes(lx, ly, 8, 4);
  const SpectralBasis basis(lx, ly, 2, 1);
  OptimizerSettings settings;
  settings.max_strain = 0.0;
  BOOST_CHECK_THROW(
    static_cast<void>(optimize(nodes, basis, settings)),
    std::invalid_argument);
  settings.max_strain = 4.0;
  settings.stall_tol = -1.0;
  BOOST_CHECK_THROW(
    static_cast<void>(optimize(nodes, basis, settings)),
    std::invalid_argument);
}

// The flow and the modal derivatives give bitwise identical results on one
// thread and on several. Without OpenMP, both runs are serial.
BOOST_AUTO_TEST_CASE(results_do_not_depend_on_thread_count)
{
  const ThreadCountGuard guard;
  BOOST_TEST_MESSAGE(
    "OpenMP " << (openmp_enabled() ? "enabled" : "disabled") << ", "
              << max_threads() << " threads by default");

  // Three blocks of the modal sums
  const NodeSet start = asymmetric_map_nodes(lx, ly, 144, 72);
  BOOST_TEST_REQUIRE(start.size() > 2 * conformality::modal_block_size);
  const SpectralBasis basis(lx, ly, 6, 3);
  std::mt19937 rng(7);
  const auto c = random_coefficients(basis.n_modes(), 2e-3, rng);

  auto run = [&](const int threads) {
    set_max_threads(threads);
    NodeSet nodes = start;
    const auto stats = advance_nodes(nodes, basis, c, 1.0, 3);
    BOOST_TEST_REQUIRE(stats.converged);
    return std::make_pair(nodes, modal_derivatives(nodes, basis));
  };
  const auto [serial_nodes, serial_f] = run(1);
  const auto [parallel_nodes, parallel_f] = run(4);

  std::size_t position_mismatches = 0;
  std::size_t jacobian_mismatches = 0;
  for (std::size_t j = 0; j < start.size(); ++j) {
    const auto &a = serial_nodes.y[j];
    const auto &b = parallel_nodes.y[j];
    position_mismatches += same_bits({a.x, a.y}, {b.x, b.y}) ? 0 : 1;
    const auto &ja = serial_nodes.jac[j];
    const auto &jb = parallel_nodes.jac[j];
    jacobian_mismatches += same_bits(
                             {ja.a11, ja.a12, ja.a21, ja.a22},
                             {jb.a11, jb.a12, jb.a21, jb.a22})
                             ? 0
                             : 1;
  }
  BOOST_TEST(position_mismatches == 0u);
  BOOST_TEST(jacobian_mismatches == 0u);
  BOOST_TEST(same_bits(serial_f, parallel_f));
}
