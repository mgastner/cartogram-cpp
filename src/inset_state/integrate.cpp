#include "inset_state.hpp"

void InsetState::preprocess()
{
  timer.start("Preprocessing");

  // Remove tiny polygons below threshold
  if (args_.remove_tiny_polygons) {
    remove_tiny_polygons(args_.min_polygon_area);
  }

  // Rescale map to fit into a rectangular box [0, lx] * [0, ly]
  rescale_map();

  // Store original coordinates for morphing animation
  if (args_.redirect_exports_to_stdout) {
    store_original_geo_divs();
  }

  if (args_.simplify) {
    std::cerr << "Start of initial simplification of " << pos_ << std::endl;

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
    initialize_identity_proj();
    write_map(inset_name_ + "_input", args_.plot_grid, true);
  }

  timer.stop("Preprocessing");
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

void InsetState::cleanup_after_integration()
{
  // Destory FFTW plans and free memory for rho and flux initializations
  destroy_fftw_plans_for_rho();
  destroy_fftw_plans_for_flux();
  ref_to_rho_init().free();
  ref_to_rho_ft().free();
  ref_to_fluxx_init().free();
  ref_to_fluxy_init().free();
}

void InsetState::integrate(ProgressTracker &progress_tracker, size_t n_insets)
{

  timer.start(inset_name_);

  // Prepare Inset for Cartogram Generation
  // -- Set up Fourier transforms
  // -- Store initial parameters
  // -- Normlize target area
  // -- Set area errors
  prepare_for_integration();
  // progress_tracker.print_progress_mid_integration(
  //   max_area_error().value,
  //   n_geo_divs(),
  //   n_finished_integrations_);

  timer.start("Integration");
  while (continue_integrating()) {

    update_file_prefix();
    timer.start(file_prefix_);

    // 1. Fill/Rasterize Density
    fill_with_density();

    // -- and blur it to facillitate integration.
    blur_density();

    // 2. Flatten Density
    if (!flatten_density()) {

      // Flatten density has failed. Increase blur width and try again
      timer.stop(file_prefix_);
      continue;
    }

    // 3. Project Polygon Points by Interpolating "Flattened" (Projected) Proxy
    // Geometry
    project();

    // 4. Update area errors and try again if necessary
    set_area_errors();
    adjust_grid();
    progress_tracker.print_progress_mid_integration(
      max_area_error().value,
      n_geo_divs(),
      n_finished_integrations_);
    increment_integration();

    timer.stop(file_prefix_);
  }
  timer.stop("Integration");

  // Update and display progress information
  std::cerr << "Finished integrating inset " << pos_ << std::endl;
  progress_tracker.update_and_print_progress_end_integration(n_geo_divs());

  // Write SVG for this inset, if requested
  if (args_.plot_polygons) {
    if (n_insets > 1) {
      write_map(inset_name() + "_output", args_.plot_grid, false);
    } else {
      write_map(inset_name() + "_cartogram", args_.plot_grid, false);
    }
  }

  // Project original map with cumulative projection
  if (args_.redirect_exports_to_stdout and !args_.qtdt_method) {
    fill_grid_diagonals(true);
    project_with_cum_proj();
  }

  // Free reserved memory
  cleanup_after_integration();
  timer.stop(inset_name_);
}