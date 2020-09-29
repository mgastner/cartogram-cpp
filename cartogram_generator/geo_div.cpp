#include "geo_div.h"

GeoDiv::GeoDiv(std::string i): id(i)
{
  return;
}

GeoDiv::~GeoDiv(void)
{
  std::cout << "GeoDiv object is being deleted" << std::endl;  
}

void GeoDiv::push_back(PolygonWH pgn_wh)
{
  polygonsWH.push_back(pgn_wh);
  return;
}
