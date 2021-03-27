#ifndef XY_POINT_H_
#define XY_POINT_H_

#include "cgal_typedef.h"

struct XY_Point {
  double x;
  double y;

  // Constructor, setting points to 0
  XY_Point() {
    x = 0;
    y = 0;
  }

  // Constructor, withe one value given
  XY_Point(Point p) {
    x = p.x();
    y = p.y();
  }

  // Constructor with two values given
  XY_Point(double xg, double yg) {
    x = xg;
    y = yg;
  }
};

#endif
