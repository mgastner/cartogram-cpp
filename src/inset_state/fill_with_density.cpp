#include "inset_state.hpp"

void InsetState::fill_with_density()
{
  if (args_.rays) {

    // Fill density using ray-shooting method
    fill_with_density_rays();

  } else {

    // Fill density with new clipping method
    // More precise, but slower
    fill_with_density_clip();
  }

  // Plot density map if requested
  if (args_.plot_density) {
    std::string file_name = file_prefix_ + "_unblurred_density.svg";
    write_density_image(file_name);
  }
}