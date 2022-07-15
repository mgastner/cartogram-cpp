#ifndef PWH_H
#define PWH_H

#include "cgal_typedef.h"

double pwh_area(const Polygon_with_holes &pwh);
bool compare_pwh(
  const Polygon_with_holes &pwh1,
  const Polygon_with_holes &pwh2);

#endif
