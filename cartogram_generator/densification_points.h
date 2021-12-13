#ifndef DENSIFICATION_POINTS_H_
#define DENSIFICATION_POINTS_H_

#include "inset_state.h"
#include "constants.h"
#include <CGAL/intersections.h>

std::vector<Point> densification_points(Point, Point,
                                        const unsigned int,
                                        const unsigned int);

bool points_almost_equal(Point, Point);

bool almost_equal(double, double);

Point rounded_point(Point);

double rounded_to_bicimal(double, unsigned int);

double rounded_to_decimal(double);

#endif
