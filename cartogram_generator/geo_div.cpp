#include "geo_div.h"

GeoDiv::GeoDiv(const std::string i) : id_(i)
{
  return;
}

const std::string GeoDiv::id() const
{
  return id_;
}

unsigned int GeoDiv::n_polygons_with_holes() const
{
  return polygons_with_holes_.size();
}

const std::vector<Polygon_with_holes> GeoDiv::polygons_with_holes() const
{
  return polygons_with_holes_;
}

std::vector<Polygon_with_holes> *GeoDiv::ref_to_polygons_with_holes()
{
  return &polygons_with_holes_;
}

void GeoDiv::push_back(const Polygon_with_holes pgn_wh)
{
  polygons_with_holes_.push_back(pgn_wh);
  return;
}

int GeoDiv::n_points()
{
  int n_points = 0;
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

void GeoDiv::adjacent_to(const std::string id)
{
  adjacent_geodivs_.insert(id);
}

const std::set<std::string> GeoDiv::adjacent_geodivs() const
{
  return adjacent_geodivs_;
}
