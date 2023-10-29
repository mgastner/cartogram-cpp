#ifndef COMPARE_INSETS_H_
#define COMPARE_INSETS_H_

#include "inset_state.h"

double frechet_distance(const Polygon_with_holes &, const Polygon_with_holes &);
double hausdorff_distance(const Polygon_with_holes &, const Polygon_with_holes &);
double symmetric_distance(const Polygon_with_holes &, const Polygon_with_holes &);

#endif
