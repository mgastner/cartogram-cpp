#include "../constants.h"
#include "round_point.h"
#include <bit>

// TODO: THE COMMENTED-OUT CODE TREATED DIFFERENCES OF 1e-13 TO BE
//       DISTINGUISHABLE, WHICH WAS TOO STRICT IN PRACTICE. AS A TEMPORARY
//       FIX, dbl_resolution IS HARDCODED. WOULD IT BE POSSIBLE TO REMOVE
//       almost_equal() IF WE USE Simple_cartesian INSTEAD OF EPICK?
// Use machine epsilon (defined in constants.h) to get almost equal doubles.
// From https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
// bool almost_equal(const double a, const double b) {
//   return std::fabs(a - b) <= dbl_epsilon * std::fabs(a + b) * 2 ||
//          std::fabs(a - b) < dbl_min;
// }
bool almost_equal(const double a, const double b) {
  return std::fabs(a - b) <= dbl_resolution;
}

// Determine whether points are almost equal
bool points_almost_equal(const Point a, const Point b) {
  return (almost_equal(a.x(), b.x()) && almost_equal(a.y(), b.y()));
}
bool xy_points_almost_equal(const XYPoint a, const XYPoint b) {
  return (almost_equal(a.x, b.x) && almost_equal(a.y, b.y));
}

// Function to round a double to a nearby bicimal. Bicimals are a more
// "native" way of rounding than decimals because doubles are represented
// as powers of 2.
double rounded_to_bicimal(const double d,
                          const unsigned int lx,
                          const unsigned int ly)
{
  double whole;
  const double fractional = std::modf(d, &whole);
  const unsigned int n_bicimals = 50 - std::bit_width(std::max(lx, ly));
  const unsigned long int power_of_2 =
    (static_cast<unsigned long int>(1) << n_bicimals);
  const double bicimals = std::round(fractional * power_of_2)  / power_of_2;
  return whole + bicimals;
}

Point rounded_point(const Point a,
                    const unsigned int lx,
                    const unsigned int ly)
{
  return Point(rounded_to_bicimal(a.x(), lx, ly),
               rounded_to_bicimal(a.y(), lx, ly));
}

XYPoint rounded_XYpoint(const XYPoint a,
                        const unsigned int lx,
                        const unsigned int ly)
{
  return XYPoint(rounded_to_bicimal(a.x, lx, ly),
                 rounded_to_bicimal(a.y, lx, ly));
}
