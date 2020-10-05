#include "geo_div.h"

GeoDiv::GeoDiv(const std::string i): id(i)
{
  return;
}

void GeoDiv::push_back(const PolygonWH pgn_wh)
{
  polygonsWH.push_back(pgn_wh);
  return;
}
