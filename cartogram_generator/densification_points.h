#ifndef DENSIFICATION_POINTS_H_
#define DENSIFICATION_POINTS_H_

#include "inset_state.h"
#include "constants.h"
#include <CGAL/intersections.h>

std::vector<Point> densification_points(const Point,
                                        const Point,
                                        const unsigned int,
                                        const unsigned int);
bool almost_equal(const double, const double);
bool points_almost_equal(const Point, const Point);
Point rounded_point(const Point, const unsigned int, const unsigned int);

#endif
