#include "geo_div.h"

GeoDiv::GeoDiv(const std::string i): id(i)
{
  return;
}

void GeoDiv::push_back(const PolygonWH pgn_wh)
{
  polygons_with_holes.push_back(pgn_wh);
  return;
}

int GeoDiv::n_polygons_with_holes(void)
{
  return polygons_with_holes.size();
}

PolygonWH GeoDiv::get_polygon_with_hole(const unsigned int i)
{
  return polygons_with_holes[i];
}
