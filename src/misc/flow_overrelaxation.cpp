#include "flow_overrelaxation.hpp"
#include "constants.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>

namespace
{

// Centered log ratios of the region areas to their targets:
// d = ln(1 + e) minus its mean over the regions
std::vector<double> centered_log_ratios(const std::vector<double> &errors)
{
  std::vector<double> d;
  if (errors.empty()) {
    return d;
  }
  d.reserve(errors.size());
  double mean = 0.0;
  for (const double e : errors) {
    d.push_back(std::log1p(e));
    mean += d.back();
  }
  mean /= static_cast<double>(errors.size());
  for (double &x : d) {
    x -= mean;
  }
  return d;
}

}  // namespace

double next_flow_overrelaxation(
  const std::vector<double> &errors_before,
  const std::vector<double> &errors_after,
  const double omega_used)
{
  assert(errors_before.size() == errors_after.size());

  // Least-squares slope of the centered log ratios after the integration on
  // those before it. Both have mean 0, so the fit needs no intercept.
  const std::vector<double> d_before = centered_log_ratios(errors_before);
  const std::vector<double> d_after = centered_log_ratios(errors_after);
  double sum_bb = 0.0;
  double sum_ba = 0.0;
  for (std::size_t i = 0; i < d_before.size(); ++i) {
    sum_bb += d_before[i] * d_before[i];
    sum_ba += d_before[i] * d_after[i];
  }
  const double t_reached = 1.0 + omega_used;
  if (!(sum_bb > 0.0) || !(t_reached > 0.0)) {
    return 0.0;
  }
  const double beta = sum_ba / sum_bb;

  // If the integration did not reduce the errors, the share left would be 1
  // or more, and over-relaxing the next integration would amplify the
  // growth of the errors
  if (!(beta < 1.0)) {
    return 0.0;
  }

  // Share that an integration ending at t = 1 would have left
  const double b =
    std::clamp((beta + omega_used) / t_reached, 0.0, max_flow_residual_share);
  return std::min(
    max_flow_overrelaxation,
    flow_overrelaxation_safety_factor * b / (1.0 - b));
}

bool area_errors_oscillate(
  const std::vector<double> &errors_before_previous,
  const std::vector<double> &errors_before,
  const std::vector<double> &errors_after)
{
  if (errors_before_previous.empty()) {
    return false;
  }
  assert(errors_before_previous.size() == errors_before.size());
  assert(errors_before.size() == errors_after.size());

  // Scalar product of the change of the errors in the integration that has
  // just finished and the change in the integration before it
  double product = 0.0;
  for (std::size_t i = 0; i < errors_before.size(); ++i) {
    product += (errors_after[i] - errors_before[i]) *
               (errors_before[i] - errors_before_previous[i]);
  }
  return product < 0.0;
}

bool max_area_error_reduction_too_small(
  const double max_error_before,
  const double max_error_after,
  const double t_reached)
{
  return max_error_after >=
         (1.0 - t_reached * min_max_area_error_reduction) * max_error_before;
}

bool flow_overrelaxation_allowed(const std::vector<double> &errors)
{
  double max_error = 0.0;
  for (const double e : errors) {
    max_error = std::max(max_error, std::abs(e));
  }
  return max_error < max_area_error_for_flow_overrelaxation;
}
