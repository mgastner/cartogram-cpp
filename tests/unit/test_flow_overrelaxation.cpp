#define BOOST_TEST_MODULE test_flow_overrelaxation
#include "constants.hpp"
#include "flow_overrelaxation.hpp"
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <vector>

namespace tt = boost::test_tools;

namespace
{

const std::vector<double> errors_before{-0.5, 0.2, 1.5, -0.05};

std::vector<double> scaled(const double factor)
{
  std::vector<double> v;
  for (const double e : errors_before) {
    v.push_back(factor * e);
  }
  return v;
}

// Log ratios ln(1 + e) of the areas to their targets before an integration
const std::vector<double> log_ratios_before{-0.5, 0.2, 1.5, -0.05};

// Errors whose log ratios are factor times log_ratios_before, plus a common
// shift. The shift changes only the normalization of the targets, which the
// centering of the log ratios removes.
std::vector<double> errors_with_log_ratios(
  const double factor,
  const double shift = 0.0)
{
  std::vector<double> v;
  for (const double x : log_ratios_before) {
    v.push_back(std::expm1(factor * x + shift));
  }
  return v;
}

// Over-relaxation that cancels the share b, reduced by the safety factor
double expected_omega(const double b)
{
  return flow_overrelaxation_safety_factor * b / (1.0 - b);
}

}  // namespace

BOOST_AUTO_TEST_CASE(plain_integration_leaving_share_b)
{
  // An integration to t = 1 that leaves 40% of every error
  BOOST_TEST(
    next_flow_overrelaxation(
      errors_with_log_ratios(1.0),
      errors_with_log_ratios(0.4),
      0.0) == expected_omega(0.4),
    tt::tolerance(1e-12));
}

BOOST_AUTO_TEST_CASE(overrelaxed_integration_is_corrected)
{
  // An integration to t = 1.25 with b = 0.4 leaves
  // 1 - 1.25 * (1 - 0.4) = 0.25 of every error
  BOOST_TEST(
    next_flow_overrelaxation(
      errors_with_log_ratios(1.0),
      errors_with_log_ratios(0.25),
      0.25) == expected_omega(0.4),
    tt::tolerance(1e-12));
}

BOOST_AUTO_TEST_CASE(integration_stopped_early_is_corrected)
{
  // An integration that stopped at t = 0.8 with b = 0.4 leaves
  // 1 - 0.8 * (1 - 0.4) = 0.52 of every error
  BOOST_TEST(
    next_flow_overrelaxation(
      errors_with_log_ratios(1.0),
      errors_with_log_ratios(0.52),
      -0.2) == expected_omega(0.4),
    tt::tolerance(1e-12));
}

BOOST_AUTO_TEST_CASE(overshoot_gives_no_overrelaxation)
{
  // The errors changed sign: the next integration ends at t = 1
  BOOST_TEST(
    next_flow_overrelaxation(
      errors_with_log_ratios(1.0),
      errors_with_log_ratios(-0.3),
      0.0) == 0.0);
}

BOOST_AUTO_TEST_CASE(large_residual_share_is_capped)
{
  // b = 0.95 is reduced to max_flow_residual_share, and the result to
  // max_flow_overrelaxation
  BOOST_TEST(
    next_flow_overrelaxation(
      errors_with_log_ratios(1.0),
      errors_with_log_ratios(0.95),
      0.0) == max_flow_overrelaxation);
}

BOOST_AUTO_TEST_CASE(zero_errors_give_no_overrelaxation)
{
  const std::vector<double> zeros(4, 0.0);
  BOOST_TEST(next_flow_overrelaxation(zeros, scaled(0.4), 0.0) == 0.0);
}

BOOST_AUTO_TEST_CASE(slope_is_least_squares_fit_of_centered_log_ratios)
{
  // Log ratios 1 and -1 before the integration, and 0.6 and 0 after it.
  // Centered, they are 1 and -1, and 0.3 and -0.3, so the slope is
  // sum(d_before * d_after) / sum(d_before^2) = 0.3
  const std::vector<double> before{std::expm1(1.0), std::expm1(-1.0)};
  const std::vector<double> after{std::expm1(0.6), 0.0};
  BOOST_TEST(
    next_flow_overrelaxation(before, after, 0.0) == expected_omega(0.3),
    tt::tolerance(1e-12));
}

BOOST_AUTO_TEST_CASE(normalization_of_targets_does_not_matter)
{
  // Scaling all areas before or after the integration by a common factor
  // leaves the estimate unchanged
  BOOST_TEST(
    next_flow_overrelaxation(
      errors_with_log_ratios(1.0, 0.3),
      errors_with_log_ratios(0.4, -0.2),
      0.0) == expected_omega(0.4),
    tt::tolerance(1e-12));
}

BOOST_AUTO_TEST_CASE(no_oscillation_without_earlier_integration)
{
  const std::vector<double> none;
  BOOST_TEST(!area_errors_oscillate(none, errors_before, scaled(0.4)));
}

BOOST_AUTO_TEST_CASE(steady_decrease_is_no_oscillation)
{
  // Every error shrinks by the same share in both integrations
  BOOST_TEST(!area_errors_oscillate(scaled(2.0), errors_before, scaled(0.5)));
}

BOOST_AUTO_TEST_CASE(overshoot_alone_is_no_oscillation)
{
  // The second integration overshoots, but the errors keep moving in the
  // same direction
  BOOST_TEST(!area_errors_oscillate(scaled(2.0), errors_before, scaled(-0.3)));
}

BOOST_AUTO_TEST_CASE(return_after_overshoot_is_oscillation)
{
  // The first integration overshoots, and the second moves the errors back
  BOOST_TEST(area_errors_oscillate(errors_before, scaled(-0.3), scaled(0.1)));
}

BOOST_AUTO_TEST_CASE(errors_alternating_between_regions_are_oscillation)
{
  // The largest error alternates between two regions, while the scalar
  // product of each region's errors before and after stays positive
  const std::vector<double> before_previous{0.10, 0.02};
  const std::vector<double> before{0.03, 0.09};
  const std::vector<double> after{0.09, 0.03};
  BOOST_TEST(area_errors_oscillate(before_previous, before, after));
}

BOOST_AUTO_TEST_CASE(unchanged_errors_are_no_oscillation)
{
  BOOST_TEST(
    !area_errors_oscillate(errors_before, errors_before, errors_before));
}

BOOST_AUTO_TEST_CASE(plain_integration_keeps_grid_criterion)
{
  // An integration that ends at t = 1 must reduce the maximum area error by
  // min_max_area_error_reduction
  const double required = 1.0 - min_max_area_error_reduction;
  BOOST_TEST(max_area_error_reduction_too_small(1.0, required, 1.0));
  BOOST_TEST(!max_area_error_reduction_too_small(1.0, required - 0.01, 1.0));
}

BOOST_AUTO_TEST_CASE(overrelaxed_integration_must_reduce_error_more)
{
  // An integration to t = 2 must reduce the maximum area error by twice the
  // fraction
  const double required = 1.0 - 2.0 * min_max_area_error_reduction;
  BOOST_TEST(max_area_error_reduction_too_small(1.0, required + 0.01, 2.0));
  BOOST_TEST(!max_area_error_reduction_too_small(1.0, required - 0.01, 2.0));
}

BOOST_AUTO_TEST_CASE(integration_stopped_early_needs_smaller_reduction)
{
  // An integration that stopped at t = 0.5 must reduce the maximum area error
  // by half the fraction
  const double required = 1.0 - 0.5 * min_max_area_error_reduction;
  BOOST_TEST(max_area_error_reduction_too_small(1.0, required + 0.01, 0.5));
  BOOST_TEST(!max_area_error_reduction_too_small(1.0, required - 0.01, 0.5));
}

BOOST_AUTO_TEST_CASE(large_area_error_prevents_overrelaxation)
{
  BOOST_TEST(!flow_overrelaxation_allowed({0.5, -0.3, 1.5}));
}

BOOST_AUTO_TEST_CASE(maximum_area_error_must_be_below_threshold)
{
  const double threshold = max_area_error_for_flow_overrelaxation;
  BOOST_TEST(!flow_overrelaxation_allowed({threshold}));
  BOOST_TEST(flow_overrelaxation_allowed({0.99 * threshold}));
}

BOOST_AUTO_TEST_CASE(region_far_below_target_allows_overrelaxation)
{
  // A region with a tenth of its target area has the area error 0.9
  BOOST_TEST(flow_overrelaxation_allowed({0.5, -0.9}));
}
