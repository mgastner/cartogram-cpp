#ifndef ELLIPSE_H_
#define ELLIPSE_H_

#include "cgal_typedef.h"

struct Ellipse {
  double semimajor;
  double semiminor;
  Point center;
  double theta;
};

#endif
