#ifndef CONFORMALITY_FLOW_HPP_
#define CONFORMALITY_FLOW_HPP_

#include "conformality/node_set.hpp"
#include "conformality/spectral_basis.hpp"
#include <cstddef>
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

struct ImplicitMidpointSettings {

  // The fixed-point iteration stops when successive iterates differ by at
  // most residual_tol * max(lx, ly)
  double residual_tol{1e-13};
  unsigned max_iterations{50};
};

// Outcome of one implicit-midpoint substep y -> y+ that solves
// y+ = y + dt v((y + y+) / 2)
struct SubstepResult {
  Vec2 y{};

  // Exact tangent map dy+/dy of the discrete step, evaluated at the
  // converged midpoint
  Mat2 tangent{};
  unsigned iterations{0};

  // |y + dt v(ybar) - y+| at the returned point
  double residual{0.0};
  bool converged{false};
};

[[nodiscard]] SubstepResult implicit_midpoint_substep(
  const SpectralBasis &basis,
  const std::vector<double> &c,
  Vec2 y,
  double dt,
  const ImplicitMidpointSettings &settings,
  SpectralBasis::Workspace &ws);

// Summary of a multi-substep integration
struct FlowStats {
  bool converged{true};
  unsigned max_iterations{0};
  double max_residual{0.0};
};

// Advances all nodes through the flow of the fixed velocity with stream
// coefficients c for total time h, using n_substeps equal implicit-midpoint
// substeps shared by all nodes. Positions are updated, and each Jacobian is
// replaced by (tangent of the step) * J. Nodes whose solve does not
// converge are still advanced; the caller must reject the trial when
// FlowStats::converged is false.
FlowStats advance_nodes(
  NodeSet &nodes,
  const SpectralBasis &basis,
  const std::vector<double> &c,
  double h,
  std::size_t n_substeps,
  const ImplicitMidpointSettings &settings = {});

// Same flow applied to passive points. If tangents is not null, it must
// have the same length as points, and each entry is replaced by
// (tangent of the step) * entry, so that starting from the identity it
// accumulates the Jacobian of the flow map.
FlowStats advance_points(
  std::vector<Vec2> &points,
  std::vector<Mat2> *tangents,
  const SpectralBasis &basis,
  const std::vector<double> &c,
  double h,
  std::size_t n_substeps,
  const ImplicitMidpointSettings &settings = {});

}  // namespace conformality

#endif  // CONFORMALITY_FLOW_HPP_
