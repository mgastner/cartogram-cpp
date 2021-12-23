#include "geo_div.h"
#include "cgal_typedef.h"
#include "constants.h"
#include <cmath>

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

double ring_lonlat_area(const Polygon ring)
{
  double a = 0.0;
  for (unsigned long i = 1; i < ring.size() - 1; ++i) {
      double to_radians = pi / 180.0;
      Point p1 = Point(ring[i-1].x() * to_radians,
                       ring[i-1].y() * to_radians);
      Point p2 = Point(ring[i].x() * to_radians,
                       ring[i].y() * to_radians);
      Point p3 = Point(ring[i+1].x() * to_radians,
                        ring[i+1].y() * to_radians);
      
      
      // a += (p2.x() - p1.x()) * (2 + sin(p1.y()) + sin(p2.y()));
      a += (p3.x() - p1.x()) * sin(p2.y());
    }
  a = a * earth_radius * earth_radius / 2.0;
  return fabs(a);
}

double GeoDiv::area_longlat() const
{
  double area_gd = 0.0;
  for (auto pwh : polygons_with_holes()) {
    Polygon ext_ring = pwh.outer_boundary();
    area_gd += ring_lonlat_area(ext_ring);
    for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
      Polygon hole = *hci;
      area_gd -= ring_lonlat_area(hole);
    }
  }
  return area_gd;
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

void GeoDiv::push_back(const Polygon_with_holes pgn_wh)
{
  polygons_with_holes_.push_back(pgn_wh);
  return;
}

std::vector<Polygon_with_holes> *GeoDiv::ref_to_polygons_with_holes()
{
  return &polygons_with_holes_;
}
