#include "flow_overrelaxation.hpp"
#include "constants.hpp"
#include <algorithm>
#include <cassert>
#include <cstddef>

double next_flow_overrelaxation(
  const std::vector<double> &errors_before,
  const std::vector<double> &errors_after,
  const double omega_used)
{
  assert(errors_before.size() == errors_after.size());

  // Least-squares slope of the errors after the integration on the errors
  // before it, through the origin
  double sum_bb = 0.0;
  double sum_ba = 0.0;
  for (std::size_t i = 0; i < errors_before.size(); ++i) {
    sum_bb += errors_before[i] * errors_before[i];
    sum_ba += errors_before[i] * errors_after[i];
  }
  const double t_reached = 1.0 + omega_used;
  if (!(sum_bb > 0.0) || !(t_reached > 0.0)) {
    return 0.0;
  }
  const double beta = sum_ba / sum_bb;

  // Share that an integration ending at t = 1 would have left
  const double b =
    std::clamp((beta + omega_used) / t_reached, 0.0, max_flow_residual_share);
  return std::min(
    max_flow_overrelaxation,
    flow_overrelaxation_safety_factor * b / (1.0 - b));
}
