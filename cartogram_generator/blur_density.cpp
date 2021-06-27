#include "constants.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include "write_eps.h"
#include <iostream>

void blur_density(const double blur_width,
                  InsetState *inset_state,
                  bool trigger_write_density_to_eps)
{
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  FTReal2d &rho_ft = *inset_state->ref_to_rho_ft();
  const double prefactor = -0.5 * blur_width * blur_width * pi * pi;
  for (unsigned int i=0; i<lx; ++i) {
    const double scaled_i_squared = ((double) i / lx) * ((double) i / lx);
    for (unsigned int j=0; j<ly; ++j) {
      const double scaled_j_squared = ((double) j / ly) * ((double) j / ly);
      rho_ft(i, j) *=
        exp(prefactor * (scaled_i_squared + scaled_j_squared)) / (4*lx*ly);
    }
  }
  inset_state->execute_fftw_bwd_plan();
  if (trigger_write_density_to_eps) {
    std::string file_name =
      inset_state->inset_name() +
      "_blurred_density_" +
      std::to_string(inset_state->n_finished_integrations()) +
      ".eps";
    std::cout << "Writing " << file_name << std::endl;
    write_density_to_eps(file_name,
                         inset_state->rho_init().as_1d_array(),
                         inset_state);
  }
  return;
}
