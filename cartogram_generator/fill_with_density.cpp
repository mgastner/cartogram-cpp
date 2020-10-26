#include "map_state.h"

void fill_with_density(MapState *map_state)
{
  std::cout << "In fill_with_density()" << std::endl;
  FTReal2d &rho_init = *map_state->ref_to_rho_init();
  for (unsigned int i=0; i<map_state->lx(); i++) {
    for (unsigned int j=0; j<map_state->ly(); j++) {
      if (i >= 0.25 * map_state->lx() && i < 0.75 * map_state->lx() &&
          j >= 0.25 * map_state->ly() && j < 0.75 * map_state->ly()) {
        rho_init(i, j) = 1.0;
      } else {
        rho_init(i, j) = 0.0;
      }
    }
  }
  map_state->execute_fwd_plan();
  return;
}
