#ifndef INTERSECTION_H_
#define INTERSECTION_H_

#include <string>
#include "xy_point.h"

// Struct to store intersection between line parallel with
// axis, and line segment.
struct intersection {

  // Intersection coordinates
  // The x OR y coordinate, depending on which axis is the line parallel to.
  // The coordinate that does not represent the line is stored.
  double coord;
  double target_density;  // GeoDiv's target_density
  std::string geo_div_id;  // GeoDIv's ID
  bool direction;  // Does intersection enter (true) or exit (false)?

  // Overloading "<" operator, similar to above
  bool operator < (const intersection &rhs) const
  {
    return (coord < rhs.coord ||
            (coord == rhs.coord && direction < rhs.direction));
  }
};


bool ray_y_intersects(XYPoint a,
                      XYPoint b,
                      double ray_y,
                      intersection *temp,
                      double target_density,
                      double epsilon);

bool ray_x_intersects(XYPoint a,
                      XYPoint b,
                      double ray_x,
                      intersection *temp,
                      double target_density,
                      double epsilon);

#endif
