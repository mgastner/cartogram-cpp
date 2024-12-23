#include "constants.hpp"
#include "inset_state.hpp"

void InsetState::blur_density(const double blur_width, bool plot_density)
{
  const double prefactor = -0.5 * blur_width * blur_width * pi * pi;
#pragma omp parallel for default(none) shared(prefactor)
  for (unsigned int i = 0; i < lx_; ++i) {
    const double scaled_i = static_cast<double>(i) / lx_;
    const double scaled_i_squared = scaled_i * scaled_i;
    for (unsigned int j = 0; j < ly_; ++j) {
      const double scaled_j = static_cast<double>(j) / ly_;
      const double scaled_j_squared = scaled_j * scaled_j;
      rho_ft_(i, j) *= exp(prefactor * (scaled_i_squared + scaled_j_squared)) /
                       (4 * lx_ * ly_);
    }
  }

  execute_fftw_bwd_plan();

  // Do not plot if the blur width is too small
  if (plot_density && blur_width > 0.1) {
    std::string file_name = inset_name_ + "_blurred_density_" +
                            std::to_string(n_finished_integrations()) + ".svg";
    write_density_image(file_name, false);
  }
}
