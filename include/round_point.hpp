#ifndef ROUND_POINT_HPP_
#define ROUND_POINT_HPP_

#include "cgal_typedef.hpp"

bool almost_equal(const double, const double);
bool almost_equal(const Point &, const Point &);
bool less_than(const double, const double);
bool less_than(const Point &, const Point &);
double rounded_to_bicimal(
  const double,
  const unsigned int,
  const unsigned int);
Point rounded_point(const Point &, const unsigned int, const unsigned int);
Point rounded_point(const Point &, const unsigned int);

#endif  // ROUND_POINT_HPP_
