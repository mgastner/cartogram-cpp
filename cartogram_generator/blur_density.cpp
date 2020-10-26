#include "constants.h"
#include "map_state.h"
#include <iostream>

void blur_density(const double blur_width, MapState *map_state)
{
  std::cout << "In blur_density()" << std::endl;
  const unsigned int lx = map_state->lx();
  const unsigned int ly = map_state->ly();
  FTReal2d &rho_ft = *map_state->ref_to_rho_ft();
  const double prefactor = -0.5 * blur_width * blur_width * pi * pi;
  for (unsigned int i=0; i<lx; ++i) {
    const double scaled_i_squared = ((double) i / lx) * ((double) i / lx);
    for (unsigned int j=0; j<ly; ++j) {
      const double scaled_j_squared = ((double) j / ly) * ((double) j / ly);
      rho_ft(i, j) *=
        exp(prefactor * (scaled_i_squared + scaled_j_squared)) / (4*lx*ly);
    }
  }
  map_state->execute_bwd_plan();
  return;
}
