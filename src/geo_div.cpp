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
  for (auto pwh : polygons_with_holes()) {
    Polygon ext_ring = pwh.outer_boundary();
    a += ext_ring.area();
    for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
      Polygon hole = *hci;
      a += hole.area();
    }
  }
  return a;
}

const std::string GeoDiv::id() const
{
  return id_;
}

unsigned int GeoDiv::n_points() const
{
  unsigned int n_points = 0;
  for (Polygon_with_holes pgn_wh : polygons_with_holes_) {
    Polygon outer = pgn_wh.outer_boundary();
    n_points += outer.size();
    std::vector<Polygon> holes_v(pgn_wh.holes_begin(), pgn_wh.holes_end());
    for (Polygon hole : holes_v) {
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

const std::vector<Polygon_with_holes> GeoDiv::polygons_with_holes() const
{
  return polygons_with_holes_;
}

void GeoDiv::push_back(const Polygon_with_holes pgnwh)
{
  polygons_with_holes_.push_back(pgnwh);
  return;
}

std::vector<Polygon_with_holes> *GeoDiv::ref_to_polygons_with_holes()
{
  return &polygons_with_holes_;
}

// Function that takes a Polygon_with_holes and returns the midpoint of the
// longest line segment that is inside the polygon and is halfway through the
// northern and southern tip of the polygon. Reference:
// https://gis.stackexchange.com/questions/76498/how-is-st-pointonsurface-
// calculated
Point pgnwh_point_on_surface(const Polygon_with_holes pgnwh)
{
  // Calculate line_y
  CGAL::Bbox_2 bb = pgnwh.bbox();
  double line_y = (bb.ymin() + bb.ymax()) / 2;

  // Epsilon based on default resolution
  double epsilon = 1e-6 * (1.0/default_res);

  // Vector to store intersections
  std::vector<intersection> intersections;

  // Getting outer_boundary from pgnwh
  Polygon ext_ring = pgnwh.outer_boundary();

  // Setting up previous point to form segment (curr_point, prev_point)
  XYPoint prev_point;
  prev_point.x = ext_ring[ext_ring.size()-1][0];
  prev_point.y = ext_ring[ext_ring.size()-1][1];

  // Finding all the intersections with exterior ring
  for (unsigned int l = 0; l < ext_ring.size(); ++l) {
    XYPoint curr_point;
    curr_point.x = ext_ring[l][0];
    curr_point.y = ext_ring[l][1];
    intersection temp;

    // Function to calculate whether intersection exists between
    // segment (prev_point, curr_point) and line_y
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

  // Finding all intersections for holes
  for (auto hci = pgnwh.holes_begin(); hci != pgnwh.holes_end(); ++hci) {
    Polygon hole = *hci;
    prev_point.x = hole[hole.size()-1][0];
    prev_point.y = hole[hole.size()-1][1];
    for (unsigned int l = 0; l < hole.size(); ++l) {
      XYPoint curr_point;
      curr_point.x = hole[l][0];
      curr_point.y = hole[l][1];
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

  // Assign directions (i.e. whether the line is entering or leaving the
  // polygon with holes)
  for (unsigned int l = 0; l < intersections.size(); ++l) {
    intersections[l].direction = (l%2 == 0);
  }

  // Assigning length of line segments (to find longest) using the
  // target_density property of intersections for line segment lengths
  for (unsigned int l = 0; l < intersections.size(); l += 2) {
    intersections[l].target_density =
      intersections[l + 1].coord - intersections[l].coord;
  }

  // Finding maximum segment length
  double max_length = intersections[0].target_density;
  double left = intersections[0].coord;
  double right = intersections[1].coord;
  XYPoint midpoint((right + left) / 2, line_y);

  // Iterating through lengths
  for (unsigned int l = 0; l < intersections.size(); l += 2) { \
    if (intersections[l].target_density > max_length) {
      left = intersections[l].coord;
      right = intersections[l + 1].coord;
      max_length = intersections[l].target_density;
      midpoint.x = (right + left) / 2;
    }
  }

  // Making final midpoint and returning it
  return Point(midpoint.x, midpoint.y);
}

// Function to find point of surface of largest Polygon_with_holes
Point GeoDiv::point_on_surface() const
{

  // Finding pgnwh with largest area
  double max_area = -dbl_inf;
  Polygon_with_holes largest_pgnwh;
  for (auto pgnwh : polygons_with_holes()) {
    double area = 0.0;
    Polygon ext_ring = pgnwh.outer_boundary();
    area += ext_ring.area();
    for (auto hci = pgnwh.holes_begin(); hci != pgnwh.holes_end(); ++hci) {
      Polygon hole = *hci;
      area += hole.area();
    }
    if (area > max_area) {
      max_area = area;
      largest_pgnwh = Polygon_with_holes(pgnwh);
    }
  }

  // Returning point of surface of found polygon
  return pgnwh_point_on_surface(largest_pgnwh);
}
