#include "constants.h"
#include "inset_state.h"

#include <iostream>

void InsetState::blur_density(
  const double blur_width,
  const bool plot_density,
  const bool is_format_ps)
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

  // CSV with blurred data
  // std::string file_name =
  //   inset_name_ + "_blurred_density_" + std::to_string(n_finished_integrations());
  // std::ofstream f_csv;
  // f_csv.open(file_name + ".csv");
  // for (unsigned int i = 0; i < lx_; ++i) {
  //   for (unsigned int j = 0; j < ly_; ++j) {
  //     f_csv << i << ", " << j << ", " << rho_init_(i, j) << "\n";
  //   }
  // }

  // f_csv.close();
  execute_fftw_bwd_plan();
  if (plot_density) {
    std::string file_name = inset_name_ + "_blurred_density_" +
                            std::to_string(n_finished_integrations());
    is_format_ps ? file_name += ".ps" : file_name += ".svg";
    std::cerr << "Writing " << file_name << std::endl;
    write_density_image(
      file_name,
      rho_init_.as_1d_array(),
      false,
      is_format_ps);
  }
}
