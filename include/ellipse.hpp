#ifndef ELLIPSE_HPP_
#define ELLIPSE_HPP_

#include "cgal_typedef.hpp"

struct Ellipse {
  double semimajor;
  double semiminor;
  Point center;
  double theta;

  // To avoid calculating trigonometric functions, we do not directly store
  // the angle theta between the x-axis and the semimajor axis. Instead, we
  // store the cosine and sine of theta.
  double cos_theta;
  double sin_theta;
};

#endif  // ELLIPSE_HPP_
