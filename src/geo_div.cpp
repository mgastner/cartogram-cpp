#include "geo_div.h"

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
    const Polygon ext_ring = pwh.outer_boundary();
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
  for (const auto &pwh : polygons_with_holes_) {
    const Polygon outer = pwh.outer_boundary();
    n_points += outer.size();
    const std::vector<Polygon> holes_v(pwh.holes_begin(), pwh.holes_end());
    for (const auto &hole : holes_v) {
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

void GeoDiv::push_back(const Polygon_with_holes pwh)
{
  polygons_with_holes_.push_back(pwh);
  return;
}

std::vector<Polygon_with_holes> *GeoDiv::ref_to_polygons_with_holes()
{
  return &polygons_with_holes_;
}
