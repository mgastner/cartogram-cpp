#include "geo_div.h"

GeoDiv::GeoDiv(const std::string i): id_(i)
{
  return;
}

const int GeoDiv::n_polygons_with_holes() const
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
