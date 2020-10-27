#include "map_state.h"
#include "boost/multi_array.hpp"

void flatten_density(MapState *map_state)
{
  std::cout << "In flatten_density()" << std::endl;
  const unsigned int lx = map_state->lx();
  const unsigned int ly = map_state->ly();

  // Allocate memory for the velocity grid
  boost::multi_array<double, 2>  grid_vx(boost::extents[lx][ly]);
  boost::multi_array<double, 2>  grid_vy(boost::extents[lx][ly]);

  // Prepare Fourier transforms for the flux
  FTReal2d grid_fluxx_init(lx, ly);
  FTReal2d grid_fluxy_init(lx, ly);

  for (unsigned int i=0; i<lx; i++) {
    for (unsigned int j=0; j<ly; j++) {
      grid_fluxx_init(i, j) = 1.0;
      grid_fluxy_init(i, j) = -1.0;
    }
  }

  return;
}
