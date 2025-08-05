#include "round_point.hpp"
#include "constants.hpp"
#include <bit>

// TODO: THE COMMENTED-OUT CODE TREATED DIFFERENCES OF 1e-13 TO BE
//       DISTINGUISHABLE, WHICH WAS TOO STRICT IN PRACTICE. AS A TEMPORARY
//       FIX, dbl_resolution IS HARDCODED. WOULD IT BE POSSIBLE TO REMOVE
//       almost_equal() IF WE USE Simple_cartesian INSTEAD OF EPICK?
// Use machine epsilon (defined in constants.hpp) to get almost equal doubles.
// From https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
// bool almost_equal(const double a, const double b) {
//   return std::fabs(a - b) <= dbl_epsilon * std::fabs(a + b) * 2 ||
//          std::fabs(a - b) < dbl_min;
// }
bool almost_equal(const double a, const double b)
{
  return std::fabs(a - b) <= dbl_resolution;
}

// Determine whether points are almost equal
bool almost_equal(const Point &a, const Point &b)
{
  return (almost_equal(a.x(), b.x()) && almost_equal(a.y(), b.y()));
}

bool less_than(const double a, const double b)
{
  return !(almost_equal(a, b) || a >= b);
}

bool less_than(const Point &a, const Point &b)
{
  return !(almost_equal(a, b) || a >= b);
}
