#ifndef ROUND_POINT_HPP_
#define ROUND_POINT_HPP_

#include "cgal_typedef.hpp"

bool almost_equal(const double, const double);
bool almost_equal(const Point &, const Point &);
bool less_than(const double, const double);
bool less_than(const Point &, const Point &);

#endif  // ROUND_POINT_HPP_
