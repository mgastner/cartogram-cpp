#ifndef INTERSECTION_H_
#define INTERSECTION_H_

#include "../xy_point.h"
#include <iostream>

// Struct to store intersection between line segment and grid line.
class intersection {

 private:
  // Intersection coordinates
  // The x OR y coordinate, depending on which axis is the line parallel to.
  // The coordinate that does not represent the line is stored.
  double coord;
  bool is_x; // Is this the x coordinate (true) or the y coordinate (true)?

 public:
  double target_density;  // GeoDiv's target_density
  std::string geo_div_id;  // GeoDIv's ID
  bool direction;  // Does intersection enter (true) or exit (false)?

  // Overloading "<" operator, similar to above
  bool operator < (const intersection &rhs) const
  {
    return (coord < rhs.coord ||
            (coord == rhs.coord && direction < rhs.direction));
  }

  // Constructors
  intersection(bool);
  intersection();

  double x() const;
  double y() const;
  bool ray_intersects(XYPoint,
                      XYPoint,
                      double,
                      double,
                      double);
};

#endif
