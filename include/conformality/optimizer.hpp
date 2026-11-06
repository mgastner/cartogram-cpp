#ifndef CONFORMALITY_OPTIMIZER_HPP_
#define CONFORMALITY_OPTIMIZER_HPP_

#include "conformality/flow.hpp"
#include "conformality/node_set.hpp"
#include "conformality/spectral_basis.hpp"
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <limits>
#include <string>
#include <vector>

namespace conformality
{

// Diagnostics of one outer iteration
struct IterationRecord {
  unsigned iteration{0};

  // State at the start of the iteration
  double energy{0.0};
  double r{0.0};
  double slope{0.0};

  // Accepted trial (zero if none was accepted). The strain is h times the
  // bound on ||Dv||_F of the iteration's velocity.
  double h{0.0};
  double strain{0.0};
  std::size_t n_substeps{0};
  unsigned backtracks{0};
  double energy_after{0.0};

  // Armijo margin: E + sigma h d - E(trial), non-negative when accepted
  double armijo_margin{0.0};
  double max_det_drift{0.0};
  unsigned max_solver_iterations{0};
  double max_solver_residual{0.0};

  // Wall-clock time of the iteration in seconds
  double seconds{0.0};
};

struct OptimizerSettings {

  // Armijo constant: accept a trial of duration h if
  // E(trial) <= E + sigma h d, where d < 0 is the initial slope.
  //
  // On a quadratic model of E along the flow, Armijo admits every duration
  // up to 2 (1 - sigma) times the optimal one, and any admitted duration
  // still achieves at least 4 sigma (1 - sigma) of the optimal decrease:
  // 75% for sigma = 0.25, but only 0.04% for sigma = 1e-4. With the small
  // value and the step doubling below, the recovery tests stalled near
  // twice the optimal duration (more than 1000 iterations instead of 4).
  double sigma{0.25};

  // Duration of the first trial. If not positive, the first trial uses
  // E / R, the time at which the linear model of E reaches zero.
  double initial_step{0.0};

  // Upper bound on any trial duration
  double max_step{std::numeric_limits<double>::infinity()};

  // Upper bound on the strain h B of any trial, where B bounds ||Dv||_F
  // everywhere. It limits how much one step can deform the map (for the
  // exact flow map A_h, ||DA_h||_2 <= exp(h B)), and so caps the substep
  // count of a trial at max_strain / substep_factor. On the first real maps
  // (USA and China, 16^2 and 32^2 modes), no accepted step exceeded h B = 3.2,
  // while the first trial E / R overshot up to h B = 25 and cost most of the
  // run time.
  double max_strain{4.0};

  // Each rejected trial halves the duration, at most this many times per
  // outer iteration
  unsigned max_backtracks{30};

  unsigned max_iterations{200};

  // Substeps are chosen so that dt * (bound on ||Dv||_F) <= substep_factor
  double substep_factor{0.1};

  // If a trial would need more substeps, its duration is shortened
  std::size_t max_substeps{10000};

  // Converged when sqrt(R_k / R_0) <= gradient_tol
  double gradient_tol{1e-6};

  // Converged when an accepted step decreases E by at most
  // energy_tol * E (stagnation)
  double energy_tol{1e-12};

  // Converged when the last stall_window accepted steps together decreased
  // E by at most stall_tol * E (slow tail). A window of 0 disables the
  // test. On the first real maps, E was then within 1e-5 of its value
  // after 200 iterations.
  unsigned stall_window{10};
  double stall_tol{1e-5};

  // A trial is rejected if max_j |det J_j / g0_j - 1| exceeds this value
  double det_drift_tol{1e-9};

  ImplicitMidpointSettings solver{};

  // Called after every outer iteration with its record (e.g. to print
  // progress)
  std::function<void(const IterationRecord &)> on_iteration{};
};

enum class Outcome { converged, numerical_failure, budget_exhausted };

[[nodiscard]] std::string to_string(Outcome outcome);
std::ostream &operator<<(std::ostream &os, Outcome outcome);

// One accepted correction: the flow of the velocity with stream
// coefficients `coeffs` for duration h in n_substeps equal substeps.
// Replaying all accepted steps in order applies the same correction A to
// any point set.
struct AcceptedStep {
  std::vector<double> coeffs;
  double h{0.0};
  std::size_t n_substeps{0};
};

struct OptimizationResult {
  Outcome outcome{Outcome::budget_exhausted};
  std::string reason;
  double initial_energy{0.0};
  double final_energy{0.0};
  double initial_r{0.0};
  double final_r{0.0};
  std::vector<IterationRecord> history;
  std::vector<AcceptedStep> steps;
};

// Largest |det J_j / g0_j - 1| over the nodes
[[nodiscard]] double max_det_drift(const NodeSet &nodes);

// Minimizes E0 over maps A o T reachable by successive spectral Hamiltonian
// flows, starting from the current nodes, by H1 steepest descent with
// Armijo backtracking. On return, `nodes` holds the last accepted state.
[[nodiscard]] OptimizationResult optimize(
  NodeSet &nodes,
  const SpectralBasis &basis,
  const OptimizerSettings &settings = {});

// Applies the accepted steps of an optimization to nodes (positions and
// Jacobians) or to passive points (and optionally their tangents), in order
FlowStats replay(
  NodeSet &nodes,
  const SpectralBasis &basis,
  const std::vector<AcceptedStep> &steps,
  const ImplicitMidpointSettings &solver = {});
FlowStats replay(
  std::vector<Vec2> &points,
  std::vector<Mat2> *tangents,
  const SpectralBasis &basis,
  const std::vector<AcceptedStep> &steps,
  const ImplicitMidpointSettings &solver = {});

}  // namespace conformality

#endif  // CONFORMALITY_OPTIMIZER_HPP_
