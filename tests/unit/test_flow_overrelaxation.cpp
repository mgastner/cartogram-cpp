#define BOOST_TEST_MODULE test_flow_overrelaxation
#include "constants.hpp"
#include "flow_overrelaxation.hpp"
#include <boost/test/included/unit_test.hpp>
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
    next_flow_overrelaxation(errors_before, scaled(0.4), 0.0) ==
      expected_omega(0.4),
    tt::tolerance(1e-14));
}

BOOST_AUTO_TEST_CASE(overrelaxed_integration_is_corrected)
{
  // An integration to t = 1.25 with b = 0.4 leaves
  // 1 - 1.25 * (1 - 0.4) = 0.25 of every error
  BOOST_TEST(
    next_flow_overrelaxation(errors_before, scaled(0.25), 0.25) ==
      expected_omega(0.4),
    tt::tolerance(1e-14));
}

BOOST_AUTO_TEST_CASE(integration_stopped_early_is_corrected)
{
  // An integration that stopped at t = 0.8 with b = 0.4 leaves
  // 1 - 0.8 * (1 - 0.4) = 0.52 of every error
  BOOST_TEST(
    next_flow_overrelaxation(errors_before, scaled(0.52), -0.2) ==
      expected_omega(0.4),
    tt::tolerance(1e-14));
}

BOOST_AUTO_TEST_CASE(overshoot_gives_no_overrelaxation)
{
  // The errors changed sign: the next integration ends at t = 1
  BOOST_TEST(
    next_flow_overrelaxation(errors_before, scaled(-0.3), 0.0) == 0.0);
}

BOOST_AUTO_TEST_CASE(large_residual_share_is_capped)
{
  // b = 0.95 is reduced to max_flow_residual_share, and the result to
  // max_flow_overrelaxation
  BOOST_TEST(
    next_flow_overrelaxation(errors_before, scaled(0.95), 0.0) ==
    max_flow_overrelaxation);
}

BOOST_AUTO_TEST_CASE(zero_errors_give_no_overrelaxation)
{
  const std::vector<double> zeros(4, 0.0);
  BOOST_TEST(next_flow_overrelaxation(zeros, scaled(0.4), 0.0) == 0.0);
}

BOOST_AUTO_TEST_CASE(slope_is_least_squares_fit_through_origin)
{
  // Errors after the integration that are not proportional: the slope is
  // sum(e_before * e_after) / sum(e_before^2)
  const std::vector<double> before{1.0, -1.0};
  const std::vector<double> after{0.5, -0.1};
  const double beta = (0.5 + 0.1) / 2.0;
  BOOST_TEST(
    next_flow_overrelaxation(before, after, 0.0) == expected_omega(beta),
    tt::tolerance(1e-14));
}
