#include "cgal_typedef.h"
#include "constants.h"

bool almost_equal(double a, double b) {
  return abs(a - b) <= dbl_epsilon * abs(a + b) * 2;
}

// Overriding operator from CGAL_Point
bool Point::operator==(const Point &rhs) const
{
  return almost_equal(x(), rhs.x()) && almost_equal(y(), rhs.y());
}
