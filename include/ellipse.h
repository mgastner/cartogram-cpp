#ifndef ELLIPSE_H_
#define ELLIPSE_H_

#include "cgal_typedef.h"

struct Ellipse {
  double semimajor;
  double semiminor;
  Point center;
  double theta;

  // To avoid calculating trigonometric functions, we do not directly store
  // the angle theta between the x-axis and the semimajor axis. Instead we
  // store the cosine and sine of theta.
  double cos_theta;
  double sin_theta;
};

#endif
