#include "intersection.hpp"

intersection::intersection() = default;

intersection::intersection(bool side)
    : is_x(side), ray_enters(false)  // Temporary value
{
}

double intersection::x() const
{
  return coord;
}

double intersection::y() const
{
  return coord;
}

// TODO: THE NAME ray_intersects() SOUNDS AS IF THE FUNCTION ONLY RETURNS A
//       boolean ANSWER. HOWEVER, IT ALSO SETS coord AND target_density. IS IT
//       POSSIBLE TO MOVE THE SIDE EFFECTS INTO SEPARATE FUNCTIONS?
bool intersection::ray_intersects(
  Point a,
  Point b,
  const double ray,
  const double td,
  const double epsilon)
{
  // Flip coordinates if rays are in y-direction. The formulae below are the
  // same, except that x is replaced with y and vice versa.
  if (!is_x) {
    a = Point(a.y(), a.x());
    b = Point(b.y(), b.x());
  }

  // Check whether an intersection is present
  if (
    ((a.y() <= ray && b.y() >= ray) || (a.y() >= ray && b.y() <= ray)) &&

    // Pre-condition to ignore grazing incidence (i.e., a line segment along
    // the polygon is exactly on the test ray)
    (a.y() != b.y())) {
    if (a.y() == ray) {
      a = Point(a.x(), a.y() + epsilon);
    } else if (b.y() == ray) {
      b = Point(b.x(), b.y() + epsilon);
    }

    // Edit intersection passed by reference. coord stores the x-coordinate.
    target_density = td;
    coord = (a.x() * (b.y() - ray) + b.x() * (ray - a.y())) / (b.y() - a.y());
    return true;
  }
  return false;
}

// This function adds intersections between a ray and a polygon to
// `intersections`
void add_intersections(
  std::vector<intersection> &intersections,
  const Polygon &pgn,
  const double ray,
  const double target_density,
  const double epsilon,
  const std::string &gd_id,
  const char axis)
{
  if (axis != 'x' && axis != 'y') {
    std::cerr << "Invalid axis in add_intersections()" << std::endl;
    exit(984321);
  }
  Point prev_point = pgn[pgn.size() - 1];
  for (auto p : pgn) {
    const Point curr_point = p;
    intersection temp(axis == 'x');
    if (temp.ray_intersects(
          curr_point,
          prev_point,
          ray,
          target_density,
          epsilon)) {
      temp.geo_div_id = gd_id;
      intersections.push_back(temp);
    }
    prev_point = curr_point;
  }
}
