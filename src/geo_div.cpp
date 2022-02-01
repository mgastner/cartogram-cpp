#include "geo_div.h"
#include "intersection.h"
#include "xy_point.h"
#include "constants.h"

GeoDiv::GeoDiv(const std::string i) : id_(i)
{
  return;
}

const std::set<std::string> GeoDiv::adjacent_geodivs() const
{
  return adjacent_geodivs_;
}

void GeoDiv::adjacent_to(const std::string id)
{
  adjacent_geodivs_.insert(id);
}

double GeoDiv::area() const
{
  double a = 0.0;
  for (const auto &pwh : polygons_with_holes()) {
    const auto ext_ring = pwh.outer_boundary();
    a += ext_ring.area();
    for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
      Polygon hole = *h;
      a += hole.area();
    }
  }
  return a;
}

const std::string GeoDiv::id() const
{
  return id_;
}

const Polygon_with_holes GeoDiv::largest_polygon_with_holes() const
{
  double max_area = -dbl_inf;
  Polygon_with_holes largest_pwh;
  for (const auto &pwh : polygons_with_holes_) {
    double area = 0.0;
    const auto ext_ring = pwh.outer_boundary();
    area += ext_ring.area();
    for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
      const auto hole = *h;
      area += hole.area();
    }
    if (area > max_area) {
      max_area = area;
      largest_pwh = Polygon_with_holes(pwh);
    }
  }
  return largest_pwh;
}

unsigned int GeoDiv::n_points() const
{
  unsigned int n_points = 0;
  for (const auto &pwh : polygons_with_holes_) {
    const auto outer = pwh.outer_boundary();
    n_points += outer.size();
    const std::vector<Polygon> holes(pwh.holes_begin(), pwh.holes_end());
    for (const auto &hole : holes) {
      n_points += hole.size();
    }
  }
  return n_points;
}

unsigned int GeoDiv::n_polygons_with_holes() const
{
  return polygons_with_holes_.size();
}

unsigned int GeoDiv::n_rings() const
{
  unsigned int n_rings = 0;
  for (const auto &pwh : polygons_with_holes_) {
    n_rings += pwh.number_of_holes() + 1;  // Add 1 for external ring
  }
  return n_rings;
}

// TODO: THIS IS NOT THE USUAL DEFINITION OF point_on_surface(). INSTEAD OF
//       THE LARGEST POLYGON WITH HOLES, THE LARGEST LINE SEGMENT IN ANY
//       POLYGON WITH HOLES IN THE MULTIPOLYGON IS CHOSEN.
Point GeoDiv::point_on_surface_of_geodiv() const
{
  return point_on_surface_of_polygon_with_holes(largest_polygon_with_holes());
}

// Function that takes a Polygon_with_holes and returns the midpoint of the
// longest line segment that is inside the polygon and is halfway through the
// northern and southern tip of the polygon. Reference:
// https://gis.stackexchange.com/questions/76498/how-is-st-pointonsurface-
// calculated
Point GeoDiv::point_on_surface_of_polygon_with_holes(
  const Polygon_with_holes pwh) const
{
  const auto bb = pwh.bbox();
  const double line_y = (bb.ymin() + bb.ymax()) / 2;
  const double epsilon = 1e-6 / default_max_n_graticule_rows_or_cols;

  // Vector to store intersections
  std::vector<intersection> intersections;
  const auto ext_ring = pwh.outer_boundary();

  // Store previous point to form segment (curr_point, prev_point)
  const auto last_ext_ring_pt = ext_ring[ext_ring.size()-1];
  XYPoint prev_point(last_ext_ring_pt.x(), last_ext_ring_pt.y());

  // Find all the intersections with exterior ring
  for (unsigned int l = 0; l < ext_ring.size(); ++l) {
    const XYPoint curr_point(ext_ring[l].x(), ext_ring[l].y());
    intersection temp;

    // Calculate whether intersection exists between segment (prev_point,
    // curr_point) and line_y
    if (ray_y_intersects(curr_point,
                         prev_point,
                         line_y,
                         &temp,
                         0,
                         epsilon)) {
      intersections.push_back(temp);
    }
    prev_point.x = curr_point.x;
    prev_point.y = curr_point.y;
  }

  // Find all intersections with holes
  for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
    Polygon hole = *h;
    prev_point.x = hole[hole.size()-1][0];
    prev_point.y = hole[hole.size()-1][1];
    for (unsigned int l = 0; l < hole.size(); ++l) {
      const XYPoint curr_point(hole[l].x(), hole[l].y());
      intersection temp;
      if (ray_y_intersects(curr_point,
                           prev_point,
                           line_y,
                           &temp,
                           0,
                           epsilon)) {
        intersections.push_back(temp);
      }
      prev_point.x = curr_point.x;
      prev_point.y = curr_point.y;
    }
  }
  std::sort(intersections.begin(), intersections.end());

  // Assign directions (i.e., whether the line is entering or leaving the
  // polygon with holes)
  for (unsigned int l = 0; l < intersections.size(); ++l) {
    intersections[l].direction = (l%2 == 0);
  }

  // Assign length of line segments using the target_density property of
  // intersections for line segment lengths
  for (unsigned int l = 0; l < intersections.size(); l += 2) {
    intersections[l].target_density =
      intersections[l + 1].coord - intersections[l].coord;
  }

  // Find maximum segment length
  double max_length = intersections[0].target_density;
  double left = intersections[0].coord;
  double right = intersections[1].coord;
  XYPoint midpoint((right + left) / 2, line_y);

  // Iterate over lengths
  for (unsigned int l = 0; l < intersections.size(); l += 2) { \
    if (intersections[l].target_density > max_length) {
      left = intersections[l].coord;
      right = intersections[l + 1].coord;
      max_length = intersections[l].target_density;
      midpoint.x = (right + left) / 2;
    }
  }
  return Point(midpoint.x, midpoint.y);
}

const std::vector<Polygon_with_holes> GeoDiv::polygons_with_holes() const
{
  return polygons_with_holes_;
}

void GeoDiv::push_back(const Polygon_with_holes pwh)
{
  polygons_with_holes_.push_back(pwh);
  return;
}

std::vector<Polygon_with_holes> *GeoDiv::ref_to_polygons_with_holes()
{
  return &polygons_with_holes_;
}
