#include "intersection.h"

intersection::intersection(bool side) :
  is_x(side),
  ray_enters(false)  // Temporary value
{
  return;
}

intersection::intersection() :
  is_x(true),     // Assume x if not explicitly declared
  ray_enters(false)  // Temporary value
{
  return;
}

double intersection::x() const {
  return coord;
}

double intersection::y() const {
  return coord;
}

bool intersection::ray_intersects(XYPoint a,
                                  XYPoint b,
                                  double ray,
                                  double td,
                                  double epsilon)
{
  // Flip coordinates if rays are in y-direction. The formulae below are the
  // same, except that x is replaced with y and vice versa.
  if (!is_x) {
    a.flip();
    b.flip();
  }

  // Check whether an intersection is present
  if (((a.y <= ray && b.y >= ray) ||
       (a.y >= ray && b.y <= ray)) &&

      // Pre-condition to ignore grazing incidence (i.e., a line segment along
      // the polygon is exactly on the test ray)
      (a.y != b.y)) {
    if (a.y == ray) {
      a.y += epsilon;
    } else if (b.y == ray) {
      b.y += epsilon;
    }

    // Edit intersection passed by reference. coord stores the x coordinate.
    target_density = td;
    coord = (a.x * (b.y - ray) + b.x * (ray - a.y)) / (b.y - a.y);
    return true;
  }
  return false;
}

// This function adds intersections between a ray and a polygon to
// `intersections`
void add_intersections(std::vector<intersection> &intersections,
                       Polygon pgn,
                       double ray,
                       double target_density,
                       double epsilon,
                       std::string gd_id,
                       bool is_x_axis)
{
  XYPoint prev_point;
  prev_point.x = pgn[pgn.size()-1].x();
  prev_point.y = pgn[pgn.size()-1].y();
  for (unsigned int l = 0; l < pgn.size(); ++l) {
    XYPoint curr_point;
    curr_point.x = pgn[l].x();
    curr_point.y = pgn[l].y();
    intersection temp(is_x_axis);
    if (temp.ray_intersects(curr_point,
                            prev_point,
                            ray,
                            target_density,
                            epsilon)) {
      temp.geo_div_id = gd_id;
      intersections.push_back(temp);
    }
    prev_point.x = curr_point.x;
    prev_point.y = curr_point.y;
  }
}
