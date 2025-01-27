#include "inset_state.hpp"

void InsetState::fill_with_density() {
  timer.start("Total")
  if (args_.rays) {

      // Fill density using ray-shooting method
      fill_with_density_rays();

  } else {

      // Fill density with new clipping method
      // More precise, but slower
      fill_with_density_clip();
  }
  timer.stop("Total")

  // Plot density map if requested
  if (args_.plot_density) {
      std::string file_name = inset_name_ + "_unblurred_density_" +
                              std::to_string(n_finished_integrations_) + ".svg";
      write_density_image(file_name, false);
  }
}