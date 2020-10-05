#include "map_state.h"
#include <iostream>

void rescale_map(int longer_lattice_length, MapState *map_state)
{
  std::cout << "In rescale_map()" << std::endl;
  std::cout << "Is it a world map? "
            << map_state->is_world_map()
            << std::endl;
  return;
}
