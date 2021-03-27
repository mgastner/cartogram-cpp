#ifndef FIND_INTERSECTIONS_AT_Y_H_
#define FIND_INTERSECTIONS_AT_Y_H_

#include "xy_point.h"
#include <string>

// Struct to store intersection data
struct intersection {
  double x;  // x-coordinate of intersection
  double target_density;  // GeoDiv's target_density
  std::string geo_div_id;
  bool direction;  // Does intersection enter or exit?

  // Overload "<" operator for this data type. Idea from
  // https://stackoverflow.com/questions/4892680/sorting-a-vector-of-structs
  bool operator < (const intersection &rhs) const
  {
    return (x < rhs.x || (x == rhs.x && direction < rhs.direction));
  }
};

bool line_y_intersects(XY_Point,
                       XY_Point,
                       double,
                       intersection*,
                       double,
                       unsigned int);

#endif
