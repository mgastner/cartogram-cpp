#include "inset_state.hpp"

void InsetState::preprocess()
{
  timer.start("Total");
  timer.start("Preprocessing");

  // Remove tiny polygons below threshold
  if (args_.remove_tiny_polygons) {
      remove_tiny_polygons(args_.min_polygon_area);
  }

  // Rescale map to fit into a rectangular box [0, lx] * [0, ly]
  rescale_map(
      args_.n_grid_rows_or_cols,
      args_.world);

  // Store original coordinates for morphing animation
  if (args_.redirect_exports_to_stdout) {
      store_original_geo_divs();
  }

  if (args_.simplify) {
      std::cerr << "Start of initial simplification of " << pos_
              << std::endl;

      // Simplification reduces the number of points used to represent the
      // GeoDivs in the inset, thereby reducing output file sizes and
      // run-times
      simplify(args_.target_points_per_inset);
      std::cerr << "End of initial simplification of " << pos_ << std::endl;
  }

  // Plot if requested
  if (args_.plot_polygons) {

    // Color if necessary
    auto_color();
    write_cairo_map(
        inset_name_ + "_input",
        args_.plot_grid);
  }

  timer.stop("Preprocessing");
  timer.stop("Total");
}

void InsetState::prepare_for_integration()
{
  // Prepare Inset for cartogram generation

  // Set up Fourier transforms
  ref_to_rho_init().allocate(lx_, ly_);
  ref_to_rho_ft().allocate(lx_, ly_);
  ref_to_fluxx_init().allocate(lx_, ly_);
  ref_to_fluxy_init().allocate(lx_, ly_);
  make_fftw_plans_for_rho();
  make_fftw_plans_for_flux();
  initialize_identity_proj();
  initialize_cum_proj();

  // Store initial inset area to calculate area drift,
  // set area errors based on this initial_area
  store_initial_area();
  set_area_errors();

  // Store initial target area to normalize inset areas
  store_initial_target_area();

  // Normalize total target area to be equal to initial area
  normalize_target_area();
}