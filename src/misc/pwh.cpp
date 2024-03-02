#include "pwh.h"

double pwh_area(const Polygon_with_holes &pwh)
{
  double a = 0.0;
  const Polygon &ext_ring = pwh.outer_boundary();
  a += ext_ring.area();
  for (const auto &h : pwh.holes()) {
    a += h.area();
  }
  return a;
}

// Compares two polygon with holes according to their areas
bool pwh_is_larger(
  const Polygon_with_holes &pwh1,
  const Polygon_with_holes &pwh2)
{
  return pwh_area(pwh1) > pwh_area(pwh2);
}
