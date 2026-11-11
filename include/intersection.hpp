#ifndef INTERSECTION_HPP_
#define INTERSECTION_HPP_

#include "cgal_typedef.hpp"
#include "round_point.hpp"

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
    return less_than(coord, rhs.coord) ||
           (almost_equal(coord, rhs.coord) && (ray_enters && !rhs.ray_enters));
  }

  // Constructors
  explicit intersection(bool side);
  intersection();

  [[nodiscard]] double x() const;
  [[nodiscard]] double y() const;
  bool ray_intersects(Point a, Point b, double ray, double td, double epsilon);
};

void add_intersections(
  std::vector<intersection> &intersections,
  const Polygon &pgn,
  double ray,
  double target_density,
  double epsilon,
  const std::string &gd_id,
  size_t pwh_idx,
  char axis);

#endif  // INTERSECTION_HPP_
