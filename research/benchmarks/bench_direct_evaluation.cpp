// D2 benchmark: cost of direct spectral evaluation in the conformality core
// at cartogram-like sizes (see README.md in this directory for how to build
// and run it, and for recorded results).
//
// Nodes: 512 x 512 cell-centre nodes of the asymmetric manufactured map on a
// 256 x 256 lattice (4 nodes per lattice cell). For each mode set, the
// program times the modal derivatives F, one evaluation sweep of v and Dv
// over all nodes, one implicit-midpoint substep of all nodes, and one outer
// optimizer iteration.

#include "conformality/flow.hpp"
#include "conformality/objective.hpp"
#include "conformality/optimizer.hpp"
#include "conformality_test_maps.hpp"
#include <chrono>
#include <cstdio>

using conformality::advance_nodes;
using conformality::FlowStats;
using conformality::h1_descent;
using conformality::modal_derivatives;
using conformality::NodeSet;
using conformality::optimize;
using conformality::OptimizerSettings;
using conformality::SpectralBasis;
using conformality::velocity_gradient_bound;
using clk = std::chrono::steady_clock;

namespace
{

double seconds(const clk::time_point a, const clk::time_point b)
{
  return std::chrono::duration<double>(b - a).count();
}

}  // namespace

int main()
{
  constexpr double side = 256.0;
  constexpr std::size_t nodes_per_side = 512;
  const NodeSet nodes = conformality_test::asymmetric_map_nodes(
    side,
    side,
    nodes_per_side,
    nodes_per_side);
  const double n = static_cast<double>(nodes.size());
  std::printf(
    "N = %zu nodes, NodeSet memory = %.1f MB (80 B per node)\n",
    nodes.size(),
    n * 80.0 / 1e6);

  double checksum = 0.0;
  for (const std::size_t m : {16u, 32u, 64u}) {
    const SpectralBasis basis(side, side, m, m);
    const double n_modes = static_cast<double>(m * m);

    const auto t0 = clk::now();
    const auto f = modal_derivatives(nodes, basis);
    const auto t1 = clk::now();

    const auto descent = h1_descent(basis, f);
    SpectralBasis::Workspace ws;
    const auto t2 = clk::now();
    for (std::size_t j = 0; j < nodes.size(); ++j) {
      checksum += basis.evaluate(descent.coeffs, nodes.y[j], ws).v.x;
    }
    const auto t3 = clk::now();

    // One substep of a trial whose duration gives dt * bound = 0.1
    const double bound = velocity_gradient_bound(basis, descent.coeffs);
    NodeSet stepped = nodes;
    const auto t4 = clk::now();
    const FlowStats stats =
      advance_nodes(stepped, basis, descent.coeffs, 0.1 / bound, 1);
    const auto t5 = clk::now();

    // One outer optimizer iteration from the same state
    NodeSet optimized = nodes;
    OptimizerSettings settings;
    settings.max_iterations = 1;
    const auto t6 = clk::now();
    const auto result = optimize(optimized, basis, settings);
    const auto t7 = clk::now();

    // Flop counts: 4 M per node for F (two rank-one updates), 7 M per node
    // for the three inner sums of one evaluation
    const double f_gflop = n * 4.0 * n_modes / 1e9;
    const double v_gflop = n * 7.0 * n_modes / 1e9;
    std::printf("\nmodes %zux%zu (M = %.0f)\n", m, m, n_modes);
    std::printf(
      "  F (modal derivatives): %7.3f s  (%5.1f GFLOP, %5.2f GFLOP/s)\n",
      seconds(t0, t1),
      f_gflop,
      f_gflop / seconds(t0, t1));
    std::printf(
      "  one v + Dv sweep:      %7.3f s  (%5.1f GFLOP, %5.2f GFLOP/s)\n",
      seconds(t2, t3),
      v_gflop,
      v_gflop / seconds(t2, t3));
    std::printf(
      "  one substep:           %7.3f s  (max %u fixed-point iterations, "
      "%s)\n",
      seconds(t4, t5),
      stats.max_iterations,
      stats.converged ? "converged" : "NOT converged");
    if (!result.history.empty()) {
      const auto &rec = result.history.front();
      std::printf(
        "  one outer iteration:   %7.3f s  (h = %.3g, %zu substeps, %u "
        "backtracks, E %.4g -> %.4g)\n",
        seconds(t6, t7),
        rec.h,
        rec.n_substeps,
        rec.backtracks,
        rec.energy,
        rec.energy_after);
    }
  }

  // Printed so that the evaluation loop cannot be optimized away
  std::printf("\nchecksum %.6e\n", checksum);
  return 0;
}
