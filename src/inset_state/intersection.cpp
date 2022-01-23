#include "intersection.h"

intersection::intersection(bool side) :
  is_x(side),
  direction(false)  // Temporary value
{
  return;
}

intersection::intersection() :
  is_x(true),     // Assuming x if not explicitly declared
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

  // Flip coordinates, in case y.
  // The formulae below would be the exact same, just with x replaced with y
  // and vice-versa.
  if (!is_x) {
    a.flip(); b.flip();
  }

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
    coord = (a.x * (b.y - ray) + b.x * (ray - a.y)) / (b.y - a.y);
    return true;
  }
  return false;
}
