#include "find_intersections_at_y.h"
#include "xy_point.h"

bool line_y_intersects(XY_Point a,
                       XY_Point b,
                       double line_y,
                       intersection *temp,
                       double target_density,
                       unsigned int resolution)
{
  double epsilon = 1e-6 * (1.0/resolution);

  // Check if intersection is present
  if (((a.y <= line_y && b.y >= line_y) ||
       (a.y >= line_y && b.y <= line_y)) &&

      // Pre-condition to ignore grazing incidence (i.e. a line segment along
      // the polygon is exactly on the test line)
      (a.y != b.y)) {
    if (a.y == line_y) {
      a.y += epsilon;
    } else if (b.y == line_y) {
      b.y += epsilon;
    }

    // Edit intersection passed by reference
    temp->x = (a.x * (b.y - line_y) + b.x * (line_y - a.y)) / (b.y - a.y);
    temp->target_density = target_density;
    temp->direction = false;  // Temporary value
    return true;
  }
  return false;
}
