#include "intersection.h"

intersection::intersection(bool side) :
  x_or_y(side),
  direction(false)  // Temporary value
{
  return;
}

double intersection::x() const {
  return coord;
}

double intersection::y() const {
  return coord;
}

bool intersection::ray_intersects(XYPoint a,
                                  XYPoint b,
                                  double ray,
                                  double td,
                                  double epsilon)
{
  // Check if intersection is present
  if (((a.y <= ray && b.y >= ray) ||
       (a.y >= ray && b.y <= ray)) &&

      // Pre-condition to ignore grazing incidence (i.e., a line segment along
      // the polygon is exactly on the test ray)
      (a.y != b.y)) {
    if (a.y == ray) {
      a.y += epsilon;
    } else if (b.y == ray) {
      b.y += epsilon;
    }

    // Edit intersection passed by reference
    // coord stores the x coordinate.
    target_density = td;
    if (x_or_y) {
      coord = (a.y * (b.x - ray) + b.y * (ray - a.x)) / (b.x - a.x);
    } else {
      coord = (a.x * (b.y - ray) + b.x * (ray - a.y)) / (b.y - a.y);
    }
    return true;
  }
  return false;
}
