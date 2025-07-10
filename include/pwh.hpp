#ifndef PWH_H
#define PWH_H

#include "cgal_typedef.hpp"

double pwh_area(const Polygon_with_holes &pwh);
bool pwh_is_larger(
  const Polygon_with_holes &pwh1,
  const Polygon_with_holes &pwh2);

#endif  // PWH_H
