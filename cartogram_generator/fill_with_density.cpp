#include "map_state.h"

void fill_with_density(MapState *map_state)
{
  std::cout << "In fill_with_density()" << std::endl;
  double *rho_init = map_state->get_rho_init();
  int lx = map_state->get_lx();
  int ly = map_state->get_ly();
  for (unsigned int i=0; i < lx; ++i) {
    for (unsigned int j=0; j < ly; ++j) {
      rho_init[i*ly + j] = i + j;
    }
  }
  fftw_execute(map_state->get_plan_fwd());

  double *rho_ft = map_state->get_rho_ft();
  for (unsigned int i=0; i < lx; ++i) {
    for (unsigned int j=0; j < ly; ++j) {
      std::cout << rho_ft[i*ly + j] << std::endl;
    }
  }
  return;
}
