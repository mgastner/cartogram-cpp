#ifndef FLOW_OVERRELAXATION_HPP_
#define FLOW_OVERRELAXATION_HPP_

#include <vector>

// Over-relaxation of the flow-based method.
//
// Each integration of the flow-based method moves the map with the velocity
// v = J / rho(t), where rho(t) = (1 - t) rho_init + t rho_bar. At t = 1, the
// blurred density is equalized, but the area errors are not: because of the
// blur and the finite grid, an integration that ends at t = 1 leaves a share
// b of each region's signed area error e = area / target - 1. Integrating to
// t_end = 1 + omega instead changes each error by about
// (1 + omega) (1 - b) times its value, so omega = b / (1 - b) would leave no
// error.
//
// After each integration, we estimate b from that integration alone, using
// the errors of all regions before and after it, and choose the
// over-relaxation omega of the next integration. The function returns
//   omega = min(max_flow_overrelaxation,
//               flow_overrelaxation_safety_factor * b / (1 - b)),
// where
//   beta = sum(e_before * e_after) / sum(e_before^2)
// is the least-squares slope through the origin of the errors after the
// integration on the errors before it (one point per region), and
//   b = (beta + omega_used) / (1 + omega_used),
// clamped to [0, max_flow_residual_share], removes the effect of the
// over-relaxation omega_used = t_reached - 1 of the integration that has just
// finished (negative if it stopped before t = 1). The safety factor is below
// 1 because the next integration can leave a different share than the
// integration that has just finished.
double next_flow_overrelaxation(
  const std::vector<double> &errors_before,
  const std::vector<double> &errors_after,
  double omega_used);

// Whether an integration has reduced the maximum area error too little, so
// that the grid should be refined. An integration that ends at t = 1 must
// reduce the maximum area error by at least the fraction
// min_max_area_error_reduction. An over-relaxed integration reduces each
// error by about (1 + omega) times as much as an integration that ends at
// t = 1, so an integration that reached t_reached must reduce the maximum
// area error by t_reached * min_max_area_error_reduction.
bool max_area_error_reduction_too_small(
  double max_error_before,
  double max_error_after,
  double t_reached);

#endif  // FLOW_OVERRELAXATION_HPP_
