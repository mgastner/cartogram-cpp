#include "constants.h"
#include "map_state.h"
#include "write_eps.h"
#include <iostream>

void blur_density(const double blur_width, MapState *map_state)
{
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
  if (map_state->trigger_write_density_to_eps()) {
    std::string file_name =
      std::string("blurred_density_") +
      std::to_string(map_state->n_finished_integrations()) +
      ".eps";
    std::cout << "Writing " << file_name << std::endl;
    FTReal2d &rho_init = *map_state->ref_to_rho_init();
    write_density_to_eps(file_name, rho_init.array(), map_state);
  }
  return;
}
