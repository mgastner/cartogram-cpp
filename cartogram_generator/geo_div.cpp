#include "geo_div.h"

GeoDiv::GeoDiv(const std::string i): id(i)
{
  return;
}

int GeoDiv::n_polygons_with_holes(void) const
{
  return polygons_with_holes.size();
}

std::vector<PolygonWH> GeoDiv::get_polygons_with_holes(void) const
{
  return polygons_with_holes;
}

std::vector<PolygonWH> *GeoDiv::ref_to_polygons_with_holes(void)
{
  return &polygons_with_holes;
}

void GeoDiv::push_back(const PolygonWH pgn_wh)
{
  polygons_with_holes.push_back(pgn_wh);
  return;
}
