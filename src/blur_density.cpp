#include "constants.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include "write_image.h"

void blur_density(const double blur_width,
                  const bool plot_density,
                  const bool image_format_ps,
                  InsetState *inset_state)
{
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  FTReal2d &rho_ft = *inset_state->ref_to_rho_ft();
  const double prefactor = -0.5 * blur_width * blur_width * pi * pi;
  for (unsigned int i = 0; i<lx; ++i) {
    const double scaled_i = static_cast<double>(i) / lx;
    const double scaled_i_squared = scaled_i * scaled_i;
    for (unsigned int j = 0; j<ly; ++j) {
      const double scaled_j = static_cast<double>(j) / ly;
      const double scaled_j_squared = scaled_j * scaled_j;
      rho_ft(i, j) *=
        exp(prefactor * (scaled_i_squared + scaled_j_squared)) / (4*lx*ly);
    }
  }
  inset_state->execute_fftw_bwd_plan();
  if (plot_density) {
    std::string file_name =
      inset_state->inset_name() +
      "_blurred_density_" +
      std::to_string(inset_state->n_finished_integrations());
    
    image_format_ps ? file_name += ".ps" : file_name += ".svg";
    std::cerr << "Writing " << file_name << std::endl;
    FTReal2d &rho_init = *inset_state->ref_to_rho_init();
    write_density_image(file_name, rho_init.as_1d_array(), false,
                      image_format_ps, inset_state);
  }
  return;
}
