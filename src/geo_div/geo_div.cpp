#include "geo_div.h"
#include "cgal_typedef.h"
#include "constants.h"
#include <utility>

GeoDiv::GeoDiv() = default;

GeoDiv::GeoDiv(std::string i) : id_(std::move(i)) {}

std::set<std::string> GeoDiv::adjacent_geodivs() const
{
  return adjacent_geodivs_;
}

void GeoDiv::adjacent_to(const std::string &id)
{
  adjacent_geodivs_.insert(id);
}

double GeoDiv::area() const
{
  double a = 0.0;
  for (const auto &pwh : polygons_with_holes()) {
    a += pwh_area(pwh);
  }
  return a;
}

std::string GeoDiv::id() const
{
  return id_;
}

Polygon_with_holes GeoDiv::largest_polygon_with_holes() const
{
  double max_area = -dbl_inf;
  Polygon_with_holes largest_pwh;
  for (const auto &pwh : polygons_with_holes_) {
    double area = 0.0;
    const auto &ext_ring = pwh.outer_boundary();
    area += ext_ring.area();
    for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
      area += h->area();
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
    const auto &outer = pwh.outer_boundary();
    n_points += outer.size();
    for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
      n_points += h->size();
    }
  }
  return n_points;
}

unsigned int GeoDiv::n_rings() const
{
  unsigned int n_rings = 0;
  for (const auto &pwh : polygons_with_holes_) {
    n_rings += pwh.number_of_holes() + 1;  // Add 1 for external ring
  }
  return n_rings;
}

// TODO: IS THIS THE USUAL DEFINITION OF point_on_surface()? SHOULD IT NOT BE
// THE LARGEST LINE SEGMENT IN ANY POLYGON WITH HOLES IN THE GEO_DIV?
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
  const Polygon_with_holes &pwh) const
{
  const auto bb = pwh.bbox();
  const double line_y = 0.5 * (bb.ymin() + bb.ymax());
  const double epsilon = 1e-6;

  // Vector to store intersections
  std::vector<intersection> intersections;
  add_intersections(
    intersections,
    pwh.outer_boundary(),
    line_y,
    0,
    epsilon,
    id_,
    'x');

  // Store hole intersections
  for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
    add_intersections(intersections, *hci, line_y, 0, epsilon, id_, 'x');
  }
  std::sort(intersections.begin(), intersections.end());

  // Assign directions (i.e., whether the line is entering or leaving the
  // polygon with holes)
  for (unsigned int i = 0; i < intersections.size(); ++i) {
    intersections[i].ray_enters = (i % 2 == 0);
  }

  // TODO: USING target_density WHEN WE REALLY MEAN LINE LENGTH FEELS LIKE A
  // BAD HACK. SHOULD WE RENAME THE DATA MEMBER target_density TO
  // value_in_geo_div?
  // Assign length of line segments using the target_density property of
  // intersections for line segment lengths
  for (unsigned int i = 0; i < intersections.size(); i += 2) {
    intersections[i].target_density =
      intersections[i + 1].x() - intersections[i].x();
  }

  // Find midpoint in maximum segment length
  double max_length = 0.0;
  double mid_x = -1.0;  // Temporary value

  // Iterate over lengths
  for (unsigned int i = 0; i < intersections.size(); i += 2) {
    if (intersections[i].target_density > max_length) {
      const double left = intersections[i].x();
      const double right = intersections[i + 1].x();
      max_length = intersections[i].target_density;
      mid_x = (right + left) / 2;
    }
  }
  return {mid_x, line_y};
}

std::vector<Polygon_with_holes> GeoDiv::polygons_with_holes() const
{
  return polygons_with_holes_;
}

void GeoDiv::push_back(const Polygon_with_holes &pwh)
{
  polygons_with_holes_.push_back(pwh);
}

std::vector<Polygon_with_holes> *GeoDiv::ref_to_polygons_with_holes()
{
  return &polygons_with_holes_;
}

void GeoDiv::sort_pwh()
{
  std::sort(
    polygons_with_holes_.begin(),
    polygons_with_holes_.end(),
    compare_pwh);
}
