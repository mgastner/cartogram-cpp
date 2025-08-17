#include "constants.hpp"
#include "inset_state.hpp"
#include <cmath>
#include <vector>

void InsetState::blur_density()
{
  timer.start("Blur");

  const double bw = blur_width();
  if (bw <= 0.0) {
    timer.stop("Blur");
    return;
  }

  const double prefactor = -0.5 * bw * bw * pi * pi;

  const double inv_lx2 = 1.0 / (double(lx_) * double(lx_));
  const double inv_ly2 = 1.0 / (double(ly_) * double(ly_));

  const double norm = 1.0 / (4.0 * double(lx_) * double(ly_));

  std::vector<double> w_x(lx_), w_y(ly_);

  for (unsigned int i = 0; i < lx_; ++i) {
    const double i2 = double(i) * double(i);
    w_x[i] = std::exp(prefactor * (i2 * inv_lx2)) * norm;
  }

  for (unsigned int j = 0; j < ly_; ++j) {
    const double j2 = double(j) * double(j);
    w_y[j] = std::exp(prefactor * (j2 * inv_ly2));
  }

  for (unsigned i = 0; i < lx_; ++i) {
    const double wi = w_x[i];
    for (unsigned int j = 0; j < ly_; ++j) {
      rho_ft_(i, j) *= wi * w_y[j];
    }
  }

  execute_fftw_bwd_plan();
  timer.stop("Blur");

  // Do not plot if the blur width is too small
  if (args_.plot_density && bw > 0.1) {
    std::string file_name = file_prefix_ + "_blurred_density.svg";
    write_density_image(file_name);
  }
}
