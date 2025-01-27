#include "constants.hpp"
#include "inset_state.hpp"

void InsetState::blur_density()
{
  timer.start("Total")
  // Figure out the blur width
  const double bw = blur_width();

  // No blur left to apply
  if (bw <= 0.0) return;

  timer.start("Blur");
  const double prefactor = -0.5 * bw * bw * pi * pi;
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
  timer.stop("Blur");
  timer.stop("Total")

  // Do not plot if the blur width is too small
  if (args_.plot_density && bw > 0.1) {
    std::string file_name = inset_name_ + "_blurred_density_" +
                            std::to_string(n_finished_integrations_) + ".svg";
    write_density_image(file_name, false);
  }
}
