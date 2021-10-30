#ifndef DENSIFICATION_POINTS_H_
#define DENSIFICATION_POINTS_H_

#include "inset_state.h"
#include "constants.h"
#include <CGAL/intersections.h>

std::vector<Point> densification_points(Point, Point);

bool point_almost_equal(Point, Point);

#endif
