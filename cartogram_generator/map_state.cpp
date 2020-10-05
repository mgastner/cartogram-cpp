#include "map_state.h"

MapState::MapState(bool world_proj) : world(world_proj)
{
  return;
}

int MapState::get_n_geo_divs(void)
{
  return n_geo_divs;
}

bool MapState::is_world_map(void)
{
  return world;
}
