#include "geo_div.h"

GeoDiv::GeoDiv(const std::string i) : id_(i)
{
  return;
}

const std::string GeoDiv::id() const
{
  return id_;
}

int GeoDiv::n_polygons_with_holes() const
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

double GeoDiv::area() const
{
  double a = 0.0;
  for (auto pwh : polygons_with_holes()) {
    Polygon ext_ring = pwh.outer_boundary();
    a += ext_ring.area();
    for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
      Polygon hole = *hci;
      a += hole.area();
      std::cout << "hole.area() = " << hole.area() << std::endl;
    }
  }
  return a;
}

Point GeoDiv::point_in_geodiv() const
{
  // Find largest polygon with hole in GeoDiv
  for (auto pwh : polygons_with_holes()) {
    Polygon ext_ring = pwh.outer_boundary();
    std::cout << "pos A" << std::endl;
    double a = ext_ring.area();
    std::cout << "pos B" << std::endl;
    for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
      Polygon hole = *hci;
      std::cout << "pos C" << std::endl;
      a += hole.area();
      std::cout << "pos D" << std::endl;
    }
    std::cout << "Polygon has area " << a << "\n";
  }
  return Point(0, 0);
}
