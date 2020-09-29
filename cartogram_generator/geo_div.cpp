#include "geo_div.h"

GeoDiv::GeoDiv(std::string i)
{
  id = i;
}

void GeoDiv::push_back(PolygonWH pgn_wh)
{
  polygonsWH.push_back(pgn_wh);
  return;  
}
