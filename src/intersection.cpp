#include "intersection.h"

bool ray_y_intersects(XYPoint a,
                      XYPoint b,
                      double ray_y,
                      intersection *temp,
                      double target_density,
                      double epsilon)
{
  // Check if intersection is present
  if (((a.y <= ray_y && b.y >= ray_y) ||
       (a.y >= ray_y && b.y <= ray_y)) &&

      // Pre-condition to ignore grazing incidence (i.e., a line segment along
      // the polygon is exactly on the test ray)
      (a.y != b.y)) {
    if (a.y == ray_y) {
      a.y += epsilon;
    } else if (b.y == ray_y) {
      b.y += epsilon;
    }

    // Edit intersection passed by reference
    // coord stores the x coordinate.
    temp->coord = (a.x * (b.y - ray_y) + b.x * (ray_y - a.y)) / (b.y - a.y);
    temp->target_density = target_density;
    temp->direction = false;  // Temporary value
    return true;
  }
  return false;
}


bool ray_x_intersects(XYPoint a,
                      XYPoint b,
                      double ray_x,
                      intersection *temp,
                      double target_density,
                      double epsilon)
{
  // Check if intersection is present
  if (((a.x <= ray_x && b.x >= ray_x) ||
       (a.x >= ray_x && b.x <= ray_x)) &&

      // Pre-condition to ignore grazing incidence (i.e. a line segment along
      // the polygon is exactly on the test ray)
      (a.x != b.x)) {
    if (a.x == ray_x) {
      a.x += epsilon;
    } else if (b.x == ray_x) {
      b.x += epsilon;
    }

    // Edit intersection passed by reference
    // coord stores the y coordinate.
    temp->coord = (a.y * (b.x - ray_x) + b.y * (ray_x - a.x)) / (b.x - a.x);
    temp->target_density = target_density;
    temp->direction = false;  // Temporary value
    return true;
  }
  return false;
}
