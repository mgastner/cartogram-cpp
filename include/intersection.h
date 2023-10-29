#ifndef INTERSECTION_H_
#define INTERSECTION_H_

#include "cgal_typedef.h"
#include "xy_point.h"

// Struct to store intersection between line segment and grid line
class intersection
{

private:
  double coord{};  // Coordinate in the direction of the ray
  bool is_x{};  // Is the direction along the x-axis or y-axis?

public:
  double target_density{};  // GeoDiv's target_density
  std::string geo_div_id;  // GeoDIv's ID
  bool ray_enters{};  // Does the ray enter (true) a GeoDiv or exit (false)?

  // Overloading "<" operator, similar to above
  bool operator<(const intersection &rhs) const
  {
    return (
      coord < rhs.coord ||
      (coord == rhs.coord && ray_enters < rhs.ray_enters));
  }

  // Constructors
  explicit intersection(bool);
  intersection();

  [[nodiscard]] double x() const;
  [[nodiscard]] double y() const;
  bool ray_intersects(XYPoint, XYPoint, double, double, double);
};

void add_intersections(
  std::vector<intersection> &,
  const Polygon &,
  double,
  double,
  double,
  const std::string &,
  char);

#endif
