#include "cartogram_info.hpp"
#include "constants.hpp"
#include "parse_arguments.hpp"
#include "progress_tracker.hpp"
#include "time_tracker.hpp"

int main(const int argc, const char *argv[])
{
  TimeTracker time_tracker;

  // Start of main function time
  time_tracker.start("Total Time");

  // Struct to store command line arguments
  Arguments args;

  // Parse command-line arguments
  argparse::ArgumentParser arguments = parsed_arguments(
    argc,
    argv,
  args);

  // Initialize cart_info. It contains all the information about the cartogram
  // that needs to be handled by functions called from main().
  CartogramInfo cart_info(args);

  // Determine name of input map based on the geo_file_name and store it
  std::string map_name = cart_info.set_map_name(args.visual_file_name);

  std::string crs = "+proj=longlat";
  // Read geometry. If the GeoJSON does not explicitly contain a "crs" field,
  // we assume that the coordinates are in longitude and latitude.
  cart_info.read_geojson(args.geo_file_name, args.make_csv, crs);

  // Project to equal area, if necessary
  // Write Output and EXIT if --output_equal_area_map is set to true
  cart_info.project_to_equal_area();

  std::cerr << "Coordinate reference system: " << crs << std::endl;
  if (arguments.is_used("visual_variable_file")) {

    // Read visual variables (e.g., area and color) from CSV
    cart_info.read_csv(arguments);
  }

  // Store total number of GeoDivs to monitor progress
  double total_geo_divs = cart_info.n_geo_divs();

  if (args.plot_polygons) {
    // Create copy of cart_info
    CartogramInfo tmp_ci = cart_info;

    for (InsetState &inset_state : tmp_ci.ref_to_inset_states()) {
      inset_state.normalize_inset_area(
        tmp_ci.cart_initial_total_target_area(),
        true);
      // Color if not provided
      inset_state.auto_color();
    }
    // Shift insets so that they do not overlap
    tmp_ci.reposition_insets();
    tmp_ci.write_svg("input");
  }

  // Project and exit
  if (args.output_shifted_insets) {

    // Normalize areas
    for (InsetState &inset_state : cart_info.ref_to_inset_states()) {
      if (!args.output_equal_area_map) inset_state.adjust_for_dual_hemisphere();
      inset_state.normalize_inset_area(
        cart_info.cart_initial_total_target_area(),
        true);
    }
    // Shift insets so that they do not overlap
    cart_info.reposition_insets();

    std::string suffix = "_insets_shifted";

    // Output to GeoJSON
    cart_info.write_geojson(
      args.geo_file_name,
      map_name + suffix,
      args.redirect_exports_to_stdout,
      true);
    return EXIT_SUCCESS;
  }

  // Replace missing and zero target areas with positive values
  cart_info.replace_missing_and_zero_target_areas();

  // Track progress of the cartogram generation
  ProgressTracker progress_tracker(total_geo_divs);

  // Iterate over insets
  for (InsetState &inset_state : cart_info.ref_to_inset_states()) {
    std::string inset_pos = inset_state.pos();

    // Start of inset time
    time_tracker.start("Inset " + inset_pos);

    // Determine the name of the inset
    std::string inset_name = map_name;
    if (cart_info.n_insets() > 1) {
      inset_name += "_" + inset_pos;
      std::cerr << "\nWorking on inset at position: " << inset_pos
                << std::endl;
    }
    inset_state.set_inset_name(inset_name);

    // Preprocess Inset for Integration:
    // -- Rescale
    // -- Replace missing and zero target areas,
    // -- Simplify
    // -- Remove tiny polygons
    // -- Write input map if requested (and color polygons, if necessary)
    inset_state.preprocess();

    if (args.export_preprocessed) {
      // Output rescaled GeoJSON
      cart_info.write_geojson(
        args.geo_file_name,
        // processed = simplified + rescaled
        // and potentially projected + small polygons removed
        map_name + "_input_processed");

      // Output preprocessed CSV file
      cart_info.write_csv(map_name + "_input_processed");
    }

    // Prepare Inset for Cartogram Generation
    // -- Set up Fourier transforms
    // -- Store initial parameters
    // -- Normlize target area
    // -- Set area errors
    inset_state.prepare_for_integration();

    time_tracker.start("Integration Inset " + inset_pos);
    progress_tracker.print_progress_mid_integration(inset_state);

    // Start map integration
    while (inset_state.continue_integrating()) {
      // File prefix for output files for this integration
      const std::string file_prefix =
        inset_state.inset_name() + "_" +
        std::to_string(inset_state.n_finished_integrations());

      // 1. Fill/Rasterize Density
      if (args.rays) {
        // Fill density using ray-shooting method
        time_tracker.start("Fill with Density (Ray Shooting Method)");
        inset_state.fill_with_density_rays();
        time_tracker.stop("Fill with Density (Ray Shooting Method)");
      } else {
        time_tracker.start("Fill with Density (Clipping Method)");
        inset_state.fill_with_density_clip();
        time_tracker.stop("Fill with Density (Clipping Method)");
      }

      if (args.plot_density) {
        std::string file_name = inset_state.inset_name() + "_unblurred_density_" +
                                std::to_string(inset_state.n_finished_integrations()) + ".svg";
        inset_state.write_density_image(file_name, false);
      }

      const double blur_width = inset_state.blur_width();

      if (blur_width > 0.0) {
        time_tracker.start("Blur");
        inset_state.blur_density(blur_width, args.plot_density);
        time_tracker.stop("Blur");
      }

      // 2. Flatten Density
      if (args.qtdt_method) {

        // Create Delaunay triangulation based on quadtree corners and plot
        time_tracker.start("Delaunay Triangulation");
        inset_state.create_and_store_quadtree_cell_corners();
        inset_state.create_delaunay_t();
        time_tracker.stop("Delaunay Triangulation");
        if (args.plot_quadtree) {
          inset_state.write_quadtree(file_prefix + "_quadtree");
          inset_state.write_delaunay_triangles(
            file_prefix + "a_delaunay_t",
            false);
        }

        time_tracker.start("Flatten Density (Quadtree Method)");
        if (!inset_state.flatten_density_with_node_vertices()) {

          // Flatten density has failed. Incrrease blur width and try again
          time_tracker.stop("Flatten Density (Quadtree Method)");
          inset_state.increment_n_fails_during_flatten_density();
          continue;
        }

        // Flatten density passed.
        time_tracker.stop("Flatten Density (Quadtree Method)");
        if (args.plot_quadtree) {
          inset_state.write_delaunay_triangles(
            file_prefix + "b_delaunay_t_after_flatten",
            true);
        }
      } else {

        // Using entire grid
        time_tracker.start("Flatten Density (Full Grid Method)");
        inset_state.flatten_density();
        time_tracker.stop("Flatten Density (Full Grid Method)");
      }

      // 3. Project Polygon Points by Interpolating "Flattened" (Projected) Proxy Geometry
      if (args.qtdt_method) {

        // Update triangulation adding shorter diagonal as constraint for better shape similarity
        time_tracker.start("Update Delanuay Triangulation");
        inset_state.update_delaunay_t();
        time_tracker.stop("Update Delanuay Triangulation");

        if (args.simplify) {
          time_tracker.start("Densification (using Delanuay Triangles)");
          inset_state.densify_geo_divs_using_delaunay_t();
          time_tracker.stop("Densification (using Delanuay Triangles)");
        }

        if (args.plot_quadtree) {
          inset_state.write_delaunay_triangles(
            file_prefix + "c_updated_delaunay_t_after_flatten",
            false);
        }

        // Project using the updated Delaunay triangulation and plot
        time_tracker.start("Project (Delanuay Triangulation)");
        inset_state.project_with_delaunay_t(args.redirect_exports_to_stdout);
        time_tracker.stop("Project (Delanuay Triangulation)");
        if (args.plot_quadtree) {
          inset_state.write_delaunay_triangles(
            file_prefix + "d_projected_with_updated_delaunay_t",
            true);
        }

      } else if (args.triangulation) {

        // Only densify if we will also simplify later.
        if (args.simplify) {

          // Choose diagonals that are inside grid cells, then densify.
          time_tracker.start("Densification (using Grid Diagonals)");
          inset_state.fill_grid_diagonals();
          inset_state.densify_geo_divs();
          time_tracker.stop("Densification (using Grid Diagonals)");
        }

        // Project with triangulation
        time_tracker.start("Project (Triangulation)");
        inset_state.project_with_triangulation();
        time_tracker.stop("Project (Triangulation)");
      } else {

        // Project using bilinear interpolation
        time_tracker.start("Project (Bilinear Interpolation)");
        inset_state.project();
        time_tracker.stop("Project (Bilinear Interpolation)");
      }
      if (args.simplify) {

        time_tracker.start("Simplification");
        inset_state.simplify(args.target_points_per_inset);
        time_tracker.stop("Simplification");
      }
      if (args.plot_intersections) {
        inset_state.write_intersections_image(intersections_resolution);
      }

      // 4. Update area errors and try again if necessary
      inset_state.set_area_errors();
      inset_state.adjust_grid();
      progress_tracker.print_progress_mid_integration(inset_state);
      inset_state.increment_integration();
    }
    time_tracker.stop("Integration Inset " + inset_pos);

    // Update and display progress information
    std::cerr << "Finished inset " << inset_pos << std::endl;
    progress_tracker.update_and_print_progress_end_integration(inset_state);

    if (args.plot_polygons) {
      inset_state.write_cairo_map(
        inset_state.inset_name() + "_output",
        args.plot_grid);
    }

    if (args.redirect_exports_to_stdout and !args.qtdt_method) {
      inset_state.fill_grid_diagonals(true);
      inset_state.project_with_cum_proj();
    }

    // Clean up after finishing all Fourier transforms for this inset
    inset_state.destroy_fftw_plans_for_rho();
    inset_state.destroy_fftw_plans_for_flux();
    inset_state.ref_to_rho_init().free();
    inset_state.ref_to_rho_ft().free();
    inset_state.ref_to_fluxx_init().free();
    inset_state.ref_to_fluxy_init().free();

    // End of inset time
    time_tracker.stop("Inset " + inset_pos);
  }  // End of loop over insets

  // Iterate over insets and normalize areas
  for (InsetState &inset_state : cart_info.ref_to_inset_states()) {
    std::string inset_pos = inset_state.pos();

    // Rescale insets in correct proportion to each other
    inset_state.normalize_inset_area(
      cart_info.cart_initial_total_target_area(),
      false,
      false);
    inset_state.normalize_inset_area(
      cart_info.cart_initial_total_target_area(),
      false,
      args.redirect_exports_to_stdout);
  }

  // Shift insets so that they do not overlap
  cart_info.reposition_insets(args.redirect_exports_to_stdout);

  // Output to GeoJSON
  cart_info.write_geojson(
    args.geo_file_name,
    map_name + "_cartogram",
    args.redirect_exports_to_stdout);

  if (args.plot_polygons) cart_info.write_svg("cartogram");

  // Stop of main function time
  time_tracker.stop("Total Time");

  // Print summary report
  time_tracker.print_summary_report();

  return EXIT_SUCCESS;
}