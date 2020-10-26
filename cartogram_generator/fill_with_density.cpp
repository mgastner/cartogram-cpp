#include "map_state.h"

void fill_with_density(MapState *map_state)
{
  std::cout << "In fill_with_density()" << std::endl;
  FTReal2d &rho_init = *map_state->ref_to_rho_init();
  for (unsigned int i=0; i<map_state->lx(); i++) {
    for (unsigned int j=0; j<map_state->ly(); j++) {
      rho_init(i, j) = i + j;
    }
  }
  map_state->execute_fwd_plan();

  const FTReal2d &rho_ft = *map_state->ref_to_rho_ft();
  std::cout << "Fourier transform:" << std::endl;
  for (unsigned int i=0; i<map_state->lx(); i++) {
    for (unsigned int j=0; j<map_state->ly(); j++) {
      std::cout << rho_ft(i, j) << " ";
    }
    std::cout << std::endl;
  }
  map_state->execute_bwd_plan();

  return;
}
