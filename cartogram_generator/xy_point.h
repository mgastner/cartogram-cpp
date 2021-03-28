#ifndef XY_POINT_H_
#define XY_POINT_H_

#include "cgal_typedef.h"

struct XY_Point {
  double x;
  double y;

  // Constructor, setting points to 0
  XY_Point()
  {
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

  bool is_in_polygon_with_holes(Polygon_with_holes pwh)
  {
    Point pt(this->x, this->y);

    // If point is in one of the holes, the point is not in the polygon with
    // holes
    for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
      Polygon hole = *hci;
      if (hole.bounded_side(pt) != CGAL::ON_UNBOUNDED_SIDE) {
        return false;
      }
    }
    Polygon ext_ring = pwh.outer_boundary();
    ext_ring.bounded_side(pt);
    return true;
  }
};

#endif
