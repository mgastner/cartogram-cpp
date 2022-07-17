#ifndef COMPARE_INSETS_H_
#define COMPARE_INSETS_H_

#include "inset_state.h"

double frechet_distance(Polygon_with_holes, Polygon_with_holes);
double hausdorff_distance(Polygon_with_holes, Polygon_with_holes);
double symmetric_distance(Polygon_with_holes, Polygon_with_holes);

#endif