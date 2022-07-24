#include "pwh.h"

double pwh_area(const Polygon_with_holes &pwh)
{
  double a = 0.0;
  const Polygon &ext_ring = pwh.outer_boundary();
  a += ext_ring.area();
  for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
    a += h->area();
  }
  return a;
}

// Compares two polygon with holes according to their areas
bool compare_pwh(
  const Polygon_with_holes &pwh1,
  const Polygon_with_holes &pwh2)
{
  return pwh_area(pwh1) > pwh_area(pwh2);
}
