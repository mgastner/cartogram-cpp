#ifndef ROUND_POINT_H_
#define ROUND_POINT_H_

#include "cgal_typedef.h"
#include "xy_point.h"

bool almost_equal(const double, const double);
bool points_almost_equal(const Point, const Point);
bool point_less_than(const Point, const Point);
bool xy_points_almost_equal(const XYPoint, const XYPoint);
double rounded_to_bicimal(
  const double,
  const unsigned int,
  const unsigned int);
Point rounded_point(const Point, const unsigned int, const unsigned int);

#endif
