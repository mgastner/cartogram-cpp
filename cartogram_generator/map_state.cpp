#include "map_state.h"

MapState::MapState(const bool world_proj) : world(world_proj)
{
  return;
}

void MapState::push_back(const GeoDiv gd)
{
  geo_divs.push_back(gd);
  return;
}

int MapState::n_geo_divs(void)
{
  return geo_divs.size();
}

std::vector<GeoDiv> MapState::get_geo_divs(void)
{
  return geo_divs;
}

bool MapState::is_world_map(void)
{
  return world;
}
