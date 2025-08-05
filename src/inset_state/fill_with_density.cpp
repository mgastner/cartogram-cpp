#include "inset_state.hpp"

void InsetState::fill_with_density()
{
  fill_with_density_clip();

  // Plot density map if requested
  if (args_.plot_density) {
    std::string file_name = file_prefix_ + "_unblurred_density.svg";
    write_density_image(file_name);
  }
}