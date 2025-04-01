#ifndef INTERSECTION_HPP_
#define INTERSECTION_HPP_

#include "cgal_typedef.hpp"

// Struct to store intersection between line segment and grid line
class intersection
{

private:
  double coord{};  // Coordinate in the direction of the ray
  bool is_x{};  // Is the direction along the x-axis or y-axis?

public:
  double target_density{};  // GeoDiv's target_density
  std::string geo_div_id;  // GeoDIv's ID
  size_t pwh_idx;  // Polygon ID
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
  bool ray_intersects(Point, Point, double, double, double);
};

void add_intersections(
  std::vector<intersection> &,
  const Polygon &,
  double,
  double,
  double,
  const std::string &,
  size_t,
  char);

#endif  // INTERSECTION_HPP_
