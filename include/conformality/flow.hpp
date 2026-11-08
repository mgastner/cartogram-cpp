#ifndef CONFORMALITY_FLOW_HPP_
#define CONFORMALITY_FLOW_HPP_

#include "conformality/node_set.hpp"
#include "conformality/spectral_basis.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace conformality
{

// Cayley transform (I - (dt/2) D)^(-1) (I + (dt/2) D), the tangent map of an
// implicit-midpoint step whose velocity gradient at the midpoint is D. For
// traceless D its determinant is 1 in exact arithmetic.
[[nodiscard]] Mat2 cayley(const Mat2 &d, double dt) noexcept;

// Upper bound on the Frobenius norm of the velocity gradient anywhere in
// the plane: ||Dv||_F <= sum_mn |c_mn| (kx_m^2 + ky_n^2)
[[nodiscard]] double velocity_gradient_bound(
  const SpectralBasis &basis,
  const std::vector<double> &c);

// Integrators for one substep of the flow of v = (ds/dy2, -ds/dy1). All
// three are symplectic, so in the plane every substep preserves area
// exactly, whatever its size, and each returns the exact Jacobian of its
// discrete step. They differ in cost and in how closely a substep follows
// the continuous flow.
//
//   implicit_midpoint  y+ = y + dt v((y + y+) / 2). Second order and
//                      symmetric. The equation for y+ couples both
//                      coordinates, and every fixed-point iteration
//                      evaluates the double sum over all modes.
//   symplectic_euler   with s as the Hamiltonian,
//                        y2+ = y2 - dt ds/dy1(y1, y2+),
//                        y1+ = y1 + dt ds/dy2(y1, y2+).
//                      First order. The equation for y2+ is scalar, and with
//                      y1 fixed one pass over the modes reduces every double
//                      sum to a sum over n, so a Newton iteration costs
//                      O(my) instead of O(mx my). Because s vanishes on the
//                      sides, each side maps to itself, and the substep is a
//                      one-to-one map of the rectangle onto itself as long
//                      as 1 + dt d2s/dy1dy2 > 0 everywhere, which holds when
//                      dt * velocity_gradient_bound <= 2.
//   stormer_verlet     symplectic Euler for dt / 2 as above, then its
//                      adjoint for dt / 2 (implicit in y1, see
//                      EulerVariant). Second order and symmetric; a substep
//                      costs about two symplectic Euler substeps, so at half
//                      the substep count it costs about as much as
//                      symplectic Euler.
enum class Integrator { implicit_midpoint, symplectic_euler, stormer_verlet };

// Short names "midpoint", "euler" and "verlet"
[[nodiscard]] std::string to_string(Integrator integrator);

// Inverse of to_string; throws std::invalid_argument for other names
[[nodiscard]] Integrator integrator_from_string(const std::string &name);

struct FlowSettings {
  Integrator integrator{Integrator::implicit_midpoint};

  // The implicit solves (fixed-point iteration for implicit_midpoint, Newton
  // iteration for the splitting methods) stop when successive iterates differ
  // by at most residual_tol * max(lx, ly), or fail after max_iterations
  double residual_tol{1e-13};
  unsigned max_iterations{50};
};

// Outcome of one substep y -> y+
struct SubstepResult {
  Vec2 y{};

  // Exact tangent map dy+/dy of the discrete step
  Mat2 tangent{};

  // Iterations of the implicit solve (summed over both halves for
  // stormer_verlet)
  unsigned iterations{0};

  // Residual of the implicit equation at the returned point: |y + dt v(ybar)
  // - y+| for implicit_midpoint, |p+ - p + dt dH/dq(q, p+)| of the scalar
  // equation for the splitting methods (the larger of the two halves for
  // stormer_verlet)
  double residual{0.0};

  // False if the solve did not converge, a value is not finite, or (for the
  // splitting methods) the substep is not one to one at this point
  bool converged{false};
};

[[nodiscard]] SubstepResult implicit_midpoint_substep(
  const SpectralBasis &basis,
  const std::vector<double> &c,
  Vec2 y,
  double dt,
  const FlowSettings &settings,
  SpectralBasis::Workspace &ws);

// Which coordinate a symplectic Euler substep solves for implicitly.
// implicit_in_y2 is the method of Integrator::symplectic_euler. Its adjoint,
// implicit_in_y1, is
//   y1+ = y1 + dt ds/dy2(y1+, y2),  y2+ = y2 - dt ds/dy1(y1+, y2),
// and with -dt it inverts implicit_in_y2 with dt.
enum class EulerVariant { implicit_in_y2, implicit_in_y1 };

// One symplectic Euler substep. For implicit_in_y1, the negated transposed
// coefficients are formed in the workspace on every call; advance_nodes and
// advance_points form them once per call instead.
[[nodiscard]] SubstepResult symplectic_euler_substep(
  const SpectralBasis &basis,
  const std::vector<double> &c,
  Vec2 y,
  double dt,
  const FlowSettings &settings,
  SpectralBasis::Workspace &ws,
  EulerVariant variant = EulerVariant::implicit_in_y2);

// One Stormer-Verlet substep: symplectic Euler implicit_in_y2 for dt / 2,
// then implicit_in_y1 for dt / 2
[[nodiscard]] SubstepResult stormer_verlet_substep(
  const SpectralBasis &basis,
  const std::vector<double> &c,
  Vec2 y,
  double dt,
  const FlowSettings &settings,
  SpectralBasis::Workspace &ws);

// Summary of a multi-substep integration
struct FlowStats {
  bool converged{true};
  unsigned max_iterations{0};
  double max_residual{0.0};
};

// Combines the statistics of two integrations (order independent)
void accumulate(FlowStats &total, const FlowStats &part) noexcept;

// Advances all nodes through the flow of the fixed velocity with stream
// coefficients c for total time h, using n_substeps equal substeps of the
// integrator in settings, shared by all nodes. Positions are updated, and each
// Jacobian is replaced by (tangent of the step) * J. Nodes whose solve does
// not converge are still advanced; the caller must reject the trial when
// FlowStats::converged is false. Nodes are advanced in parallel if the core
// was built with OpenMP; each node is independent of the others, so the
// result does not depend on the number of threads.
FlowStats advance_nodes(
  NodeSet &nodes,
  const SpectralBasis &basis,
  const std::vector<double> &c,
  double h,
  std::size_t n_substeps,
  const FlowSettings &settings = {});

// Same flow applied to passive points. If tangents is not null, it must
// have the same length as points, and each entry is replaced by
// (tangent of the step) * entry, so that starting from the identity it
// accumulates the Jacobian of the flow map. Parallel like advance_nodes.
FlowStats advance_points(
  std::vector<Vec2> &points,
  std::vector<Mat2> *tangents,
  const SpectralBasis &basis,
  const std::vector<double> &c,
  double h,
  std::size_t n_substeps,
  const FlowSettings &settings = {});

}  // namespace conformality

#endif  // CONFORMALITY_FLOW_HPP_
