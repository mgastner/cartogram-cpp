#include "constants.h"
#include "round_point.h"
#include <bit>

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
