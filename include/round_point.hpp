#ifndef ROUND_POINT_HPP_
#define ROUND_POINT_HPP_

#include "cgal_typedef.hpp"
#include "constants.hpp"
#include <cmath>

inline bool almost_equal(double a, double b)
{
  return std::fabs(a - b) <= dbl_resolution;
}

inline bool almost_equal(const Point &a, const Point &b)
{
  return almost_equal(a.x(), b.x()) && almost_equal(a.y(), b.y());
}

// (0: equal, -1: a < b, +1: a > b)
inline int cmp(double a, double b)
{
  if (almost_equal(a, b))
    return 0;
  return (a < b) ? -1 : 1;
}

inline int cmp(const Point &a, const Point &b)
{
  int cx = cmp(a.x(), b.x());
  if (cx != 0)
    return cx;
  return cmp(a.y(), b.y());
}

inline bool less_than(double a, double b)
{
  return cmp(a, b) < 0;
}
inline bool less_than(const Point &a, const Point &b)
{
  return cmp(a, b) < 0;
}

inline bool less_than_equal(double a, double b)
{
  return cmp(a, b) <= 0;
}
inline bool less_than_equal(const Point &a, const Point &b)
{
  return cmp(a, b) <= 0;
}

#endif  // ROUND_POINT_HPP_
