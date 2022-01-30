#include "intersection.h"

intersection::intersection(bool side) :
  is_x(side),
  ray_enters(false)  // Temporary value
{
  return;
}

intersection::intersection() :
  is_x(true),     // Assume x if not explicitly declared
  ray_enters(false)  // Temporary value
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
  // Flip coordinates if rays are in y-direction. The formulae below are the
  // same, except that x is replaced with y and vice versa.
  if (!is_x) {
    a.flip();
    b.flip();
  }

  // Check whether an intersection is present
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

    // Edit intersection passed by reference. coord stores the x coordinate.
    target_density = td;
    coord = (a.x * (b.y - ray) + b.x * (ray - a.y)) / (b.y - a.y);
    return true;
  }
  return false;
}
