#ifndef INTERSECTION_H_
#define INTERSECTION_H_

#include "xy_point.h"
#include "cgal_typedef.h"
#include <iostream>

// Struct to store intersection between line segment and grid line
class intersection {

 private:
  double coord;  // Coordinate in the direction of the ray
  bool is_x;  // Is the direction along the x-axis or y-axis?

 public:
  double target_density;  // GeoDiv's target_density
  std::string geo_div_id;  // GeoDIv's ID
  bool ray_enters;  // Does the ray enter (true) a GeoDiv or exit (false)?

  // Overloading "<" operator, similar to above
  bool operator < (const intersection &rhs) const
  {
    return (coord < rhs.coord ||
            (coord == rhs.coord && ray_enters < rhs.ray_enters));
  }

  // Constructors
  intersection(bool);
  intersection();

  double x() const;
  double y() const;
  bool ray_intersects(XYPoint,
                      XYPoint,
                      const double,
                      const double,
                      const double);
};

void add_intersections(std::vector<intersection> &,
                       const Polygon,
                       const double,
                       const double,
                       const double,
                       const std::string,
                       const char);

#endif
