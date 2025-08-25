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

  if (!args_.disable_simplification_densification) {
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
  if (args_.redirect_exports_to_stdout) {
    initialize_identity_proj();
    initialize_cum_proj();
  }

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

bool InsetState::continue_integrating() const
{

  // Calculate all the necessary information to decide whether to continue
  auto [max_area_err, worst_gd] = max_area_error();

  // A GeoDiv is still above our area error threshold
  bool area_error_above_threshold =
    max_area_err > args_.max_permitted_area_error;

  // Area expansion factor is above our threshold
  // i.e. cartogram has become too big or too small
  double area_drift = area_expansion_factor() - 1.0;
  bool area_expansion_factor_above_threshold =
    std::abs(area_drift) > max_permitted_area_drift;

  // If both the above metrics are above our threshold
  bool has_converged =
    !area_error_above_threshold && !area_expansion_factor_above_threshold;

  // Make sure to not continue endlesslely: cap at max_integrations
  bool within_integration_limit = n_finished_integrations() < max_integrations;

  // We continue if we are within the integration limit and have not converged
  bool continue_integration =
    (within_integration_limit && !has_converged) ||
    (n_finished_integrations_ < args_.min_integrations);

  // Actually hasn't converged, just reached integration limit
  if (!within_integration_limit && !has_converged) {
    converge_ = false;
    std::cerr << "ERROR: Could not converge!" << std::endl;
    if (area_error_above_threshold)
      std::cerr << "Max area error above threshold!" << std::endl;
    if (area_expansion_factor_above_threshold)
      std::cerr << "Area expansion factor above threshold!" << std::endl;
  }

  // Print control output (at end of previous integration)
  std::cerr << "Max. area err: " << max_area_err << ", GeoDiv: " << worst_gd
            << std::endl;
  std::cerr << "Current Area: " << geo_div_at_id(worst_gd).area()
            << ", Target Area: " << target_area_at(worst_gd) << std::endl;
  std::cerr << "Area drift: " << area_drift * 100.0 << "%" << std::endl;

  if (timer.total_elapsed_time_in_seconds() > args_.timeout_in_seconds) {
    std::cerr << "Warning: Timeout of " << args_.timeout_in_seconds
              << " seconds reached!" << std::endl;
    return false;
  }

  if (continue_integration) {
    // Print next integration information.
    std::cerr << "\nIntegration number " << n_finished_integrations()
              << std::endl;
    std::cerr << "Dimensions : " << lx_ << " " << ly_ << std::endl;
    std::cerr << "Number of Points: " << n_points() << std::endl;
  }
  return continue_integration;
}

void InsetState::integrate(ProgressTracker &progress_tracker)
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

    if (args_.verbose || args_.export_time_report)
      timer.start(file_prefix_);

    // 1. Fill/Rasterize Density
    fill_with_density();

    // -- and blur it to facillitate integration.
    blur_density();

    // 2. Flatten Density
    if (!flatten_density()) {

      std::cerr << "Blur width is too low. Integrator failed to flatten "
                   "density. Increasing blur width and trying again."
                << std::endl;

      increment_n_fails_during_flatten_density();

      if (args_.verbose || args_.export_time_report)
        timer.stop(file_prefix_);
      continue;
    }

    // 3. Project Polygon Points by Interpolating "Flattened" (Projected) Proxy
    // Geometry
    if (!project()) {

      std::cerr << "Triangle has flipped during triangulation. Increasing "
                   "blur width and trying again."
                << std::endl;
      increment_n_fails_during_flatten_density();

      if (args_.verbose || args_.export_time_report)
        timer.stop(file_prefix_);
      continue;
    }

    // 4. Update area errors and try again if necessary
    set_area_errors();
    adjust_grid();
    progress_tracker.print_progress_mid_integration(
      max_area_error().value,
      n_geo_divs(),
      n_finished_integrations_);
    increment_integration();

    if (args_.verbose || args_.export_time_report)
      timer.stop(file_prefix_);
  }
  timer.stop("Integration");

  // Update and display progress information
  std::cerr << "Finished integrating inset " << pos_ << std::endl;
  progress_tracker.update_and_print_progress_end_integration(n_geo_divs());

  // Write SVG for this inset, if requested
  if (args_.plot_polygons) {
    write_map(inset_name() + "_output", args_.plot_grid, false);
  }

  // Free reserved memory
  cleanup_after_integration();
  timer.stop(inset_name_);
}
