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

  std::string geo_file_name, visual_file_name;

  // Default number of grid cells along longer Cartesian coordinate axis
  unsigned int max_n_grid_rows_or_cols;

  // Target number of points to retain after simplification
  unsigned int target_points_per_inset;
  bool world;  // World maps need special projections

  // If `triangulation` is true, we apply a cartogram projection method based
  // on the triangulation of grid cells. It can eliminate intersections
  // that occur when the projected grid lines are strongly curved. Only
  // use this method if the tracer points are an FTReal2d data structure.
  bool triangulation;
  bool simplify;  // Should the polygons be simplified?

  // Other boolean values that are needed to parse the command line arguments
  bool make_csv, output_equal_area, output_to_stdout, plot_density, plot_grid,
    plot_intersections, plot_polygons, plot_quadtree, remove_tiny_polygons;

  // If the proportion of the polygon area is smaller than
  // min_polygon_area * total area, then remove polygon
  double min_polygon_area;
  bool qtdt_method;  // Use Quadtree-Delaunay triangulation

  // Parse command-line arguments
  argparse::ArgumentParser arguments = parsed_arguments(
    argc,
    argv,
    geo_file_name,
    visual_file_name,
    max_n_grid_rows_or_cols,
    target_points_per_inset,
    world,
    triangulation,
    qtdt_method,
    simplify,
    make_csv,
    output_equal_area,
    output_to_stdout,
    plot_density,
    plot_grid,
    plot_intersections,
    plot_polygons,
    remove_tiny_polygons,
    min_polygon_area,
    plot_quadtree);

  // Initialize cart_info. It contains all the information about the cartogram
  // that needs to be handled by functions called from main().
  CartogramInfo cart_info(world, visual_file_name);

  // Determine name of input map and store it
  std::string map_name = geo_file_name;
  if (map_name.find_last_of("/\\") != std::string::npos) {
    map_name = map_name.substr(map_name.find_last_of("/\\") + 1);
  }
  if (map_name.find('.') != std::string::npos) {
    map_name = map_name.substr(0, map_name.find('.'));
  }
  cart_info.set_map_name(map_name);
  if (!make_csv) {

    // Read visual variables (e.g., area and color) from CSV
    try {
      cart_info.read_csv(arguments);
    } catch (const std::system_error &e) {
      std::cerr << "ERROR reading CSV: " << e.what() << " (" << e.code() << ")"
                << std::endl;
      return EXIT_FAILURE;
    } catch (const std::runtime_error &e) {

      // If there is an error, it is probably because of an invalid CSV file
      std::cerr << "ERROR reading CSV: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Read geometry. If the GeoJSON does not explicitly contain a "crs" field,
  // we assume that the coordinates are in longitude and latitude.
  std::string crs = "+proj=longlat";
  try {
    cart_info.read_geojson(geo_file_name, make_csv, crs);
  } catch (const std::system_error &e) {
    std::cerr << "ERROR reading GeoJSON: " << e.what() << " (" << e.code()
              << ")" << std::endl;
    return EXIT_FAILURE;
  }
  std::cerr << "Coordinate reference system: " << crs << std::endl;

  // Store total number of GeoDivs to monitor progress
  double total_geo_divs = cart_info.n_geo_divs();

  // Project map and ensure that all holes are inside polygons
  for (auto &[inset_pos, inset_state] : cart_info.ref_to_inset_states()) {

    // Start of inset time
    time_tracker.start("Inset " + inset_pos);

    // Check for errors in the input topology
    try {
      inset_state.check_topology();
    } catch (const std::system_error &e) {
      std::cerr << "ERROR while checking topology: " << e.what() << " ("
                << e.code() << ")" << std::endl;
      return EXIT_FAILURE;
    }

    // Can the coordinates be interpreted as longitude and latitude?
    // TODO: The "crs" field for GeoJSON files seems to be deprecated.
    //       However, in earlier specifications, the coordinate reference
    //       system used to be written in the format specified here:
    //       https://geojson.org/geojson-spec.html#coordinate-reference-system-objects.
    //       It may be a good idea to make a list of possible entries
    //       corresponding to longitude and lattitude projection.
    //       "urn:ogc:def:crs:OGC:1.3:CRS84" is one such entry.
    const Bbox bb = inset_state.bbox();
    if (
      (bb.xmin() >= -180.0 && bb.xmax() <= 180.0) &&
      (bb.ymin() >= -90.0 && bb.ymax() <= 90.0) &&
      (crs == "+proj=longlat" || crs == "urn:ogc:def:crs:OGC:1.3:CRS84")) {

      // If yes, transform the coordinates with the Albers projection if the
      // input map is not a world map. Otherwise, use the Smyth-Craster
      // projection.
      if (world) {
        inset_state.apply_smyth_craster_projection();
      } else {
        inset_state.apply_albers_projection();
      }
    } else if (output_equal_area) {
      std::cerr << "ERROR: Input GeoJSON is not a longitude-latitude map."
                << std::endl;
      return EXIT_FAILURE;
    }

    if (simplify) {
      std::cerr << "Start of initial simplification of " << inset_pos
                << std::endl;
      time_tracker.start("Simplification");

      // Simplification reduces the number of points used to represent the
      // GeoDivs in the inset, thereby reducing output file sizes and
      // run-times
      inset_state.simplify(target_points_per_inset);

      // Update time
      time_tracker.stop("Simplification");
    }
    std::cerr << "End of initial simplification of " << inset_pos << std::endl;

    // End of inset time
    time_tracker.stop("Inset " + inset_pos);
  }

  // Replace missing and zero target areas with positive values
  cart_info.replace_missing_and_zero_target_areas();

  // Project and exit
  if (output_equal_area) {

    // Normalize areas
    for (auto &[inset_pos, inset_state] : cart_info.ref_to_inset_states()) {
      inset_state.normalize_inset_area(
        cart_info.cart_initial_total_target_area(),
        output_equal_area);
    }

    // Shift insets so that they do not overlap
    cart_info.shift_insets_to_target_position();

    // Output to GeoJSON
    cart_info.write_geojson(
      geo_file_name,
      map_name + "_equal_area.geojson",
      output_to_stdout);
    return EXIT_SUCCESS;
  }

  // Track progress of the cartogram generation
  ProgressTracker progress_tracker(total_geo_divs);

  // Iterate over insets
  for (auto &[inset_pos, inset_state] : cart_info.ref_to_inset_states()) {

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

    // Rescale map to fit into a rectangular box [0, lx] * [0, ly]
    inset_state.rescale_map(max_n_grid_rows_or_cols, cart_info.is_world_map());

    if (output_to_stdout) {

      // Store original coordinates
      inset_state.store_original_geo_divs();
    }

    // Set up Fourier transforms
    const unsigned int lx = inset_state.lx();
    const unsigned int ly = inset_state.ly();
    inset_state.ref_to_rho_init().allocate(lx, ly);
    inset_state.ref_to_rho_ft().allocate(lx, ly);
    inset_state.ref_to_fluxx_init().allocate(lx, ly);
    inset_state.ref_to_fluxy_init().allocate(lx, ly);
    inset_state.make_fftw_plans_for_rho();
    inset_state.make_fftw_plans_for_flux();
    inset_state.initialize_identity_proj();
    inset_state.initialize_cum_proj();
    inset_state.set_area_errors();

    // Store initial inset area to calculate area drift
    inset_state.store_initial_area();

    // Store initial target area to normalize inset areas
    inset_state.store_initial_target_area();

    // Normalize total target area to be equal to initial area
    inset_state.normalize_target_area();

    // Automatically color GeoDivs if no colors are provided
    if (inset_state.colors_empty()) {
      inset_state.auto_color();
    }
    if (plot_polygons) {

      // Write PNG and PS files if requested by command-line option
      std::string input_filename = inset_state.inset_name();
      if (plot_grid) {
        input_filename += "_input_grid";
      } else {
        input_filename += "_input";
      }
      std::cerr << "Writing " << input_filename << std::endl;
      inset_state.write_cairo_map(input_filename, plot_grid);
    }

    // Remove tiny polygons below threshold
    if (remove_tiny_polygons) {
      inset_state.remove_tiny_polygons(min_polygon_area);
    }

    time_tracker.start("Integration Inset " + inset_pos);

    // Start map integration
    while (inset_state.n_finished_integrations() < max_integrations &&
           (inset_state.max_area_error().value > max_permitted_area_error ||
            std::abs(inset_state.area_drift() - 1.0) > max_permitted_area_drift)) {

      std::cerr << "\nIntegration number "
                << inset_state.n_finished_integrations() << std::endl;
      std::cerr << "Number of Points: " << inset_state.n_points() << std::endl;
      if (qtdt_method) {
        time_tracker.start("Delaunay Triangulation");

        // Create the Delaunay triangulation
        inset_state.create_delaunay_t();

        time_tracker.stop("Delaunay Triangulation");

        if (plot_quadtree) {
          const std::string quadtree_filename =
            inset_state.inset_name() + "_" +
            std::to_string(inset_state.n_finished_integrations()) +
            "_quadtree";
          std::cerr << "Writing " << quadtree_filename << ".svg" << std::endl;

          // Draw the resultant quadtree
          inset_state.write_quadtree(quadtree_filename);

          const std::string delaunay_t_filename =
            inset_state.inset_name() + "_" +
            std::to_string(inset_state.n_finished_integrations()) +
            "_delaunay_t";
          std::cerr << "Writing " << delaunay_t_filename << ".svg"
                    << std::endl;
          inset_state.write_delaunay_triangles(delaunay_t_filename);
        }
      }

      const double blur_width = inset_state.blur_width();

      std::cerr << "blur_width = " << blur_width << std::endl;

      time_tracker.start("Fill with Density");

      inset_state.fill_with_density(plot_density);

      time_tracker.stop("Fill with Density");

      if (blur_width > 0.0) {
        inset_state.blur_density(blur_width, plot_density);
      }

      time_tracker.start("Flatten Density");

      if (qtdt_method) {
        inset_state.flatten_density_with_node_vertices();
      } else {
        inset_state.flatten_density();
      }

      time_tracker.stop("Flatten Density");

      if (qtdt_method) {
        if (simplify) {
          time_tracker.start("Densification");

          inset_state.densify_geo_divs_using_delaunay_t();

          time_tracker.stop("Densification");
        }

        // Project using the Delaunay triangulation
        inset_state.project_with_delaunay_t(output_to_stdout);
      } else if (triangulation) {
        time_tracker.start("Densification");

        // Choose diagonals that are inside grid cells
        inset_state.fill_grid_diagonals();

        // Densify map
        inset_state.densify_geo_divs();

        time_tracker.stop("Densification");

        // Project with triangulation
        inset_state.project_with_triangulation();
      } else {
        inset_state.project();
      }
      if (simplify) {
        time_tracker.start("Simplification");
        inset_state.simplify(target_points_per_inset);

        time_tracker.stop("Simplification");
      }
      if (plot_intersections) {
        inset_state.write_intersections_image(intersections_resolution);
      }

      // Print area drift information
      std::cerr << "Area drift: " << (inset_state.area_drift() - 1.0) * 100.0
                << "%" << std::endl;

      // Update area errors
      inset_state.set_area_errors();
      inset_state.adjust_grid();
      std::cerr << "max. area err: " << inset_state.max_area_error().value
                << ", GeoDiv: " << inset_state.max_area_error().geo_div
                << std::endl;
      progress_tracker.print_progress_mid_integration(inset_state);
      inset_state.increment_integration();
    }

    // End of inset integrations
    inset_state.check_completion(); // prints error message if conditions still not met
    time_tracker.stop("Integration Inset " + inset_pos);

    // Update and display progress information
    std::cerr << "Finished inset " << inset_pos << std::endl;
    progress_tracker.update_and_print_progress_end_integration(inset_state);

    if (plot_polygons) {
      std::string output_filename = inset_state.inset_name();
      if (plot_grid) {
        output_filename += "_output_grid";
      } else {
        output_filename += "_output";
      }
      std::cerr << "Writing " << output_filename << std::endl;
      inset_state.write_cairo_map(output_filename, plot_grid);
    }

    if (world) {
      std::string output_file_name =
        map_name + "_cartogram_in_smyth_projection.geojson";
      cart_info.write_geojson(
        geo_file_name,
        output_file_name,
        output_to_stdout);
      inset_state.revert_smyth_craster_projection();
    }

    if (output_to_stdout and !qtdt_method) {
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
  for (auto &[inset_pos, inset_state] : cart_info.ref_to_inset_states()) {

    // Rescale insets in correct proportion to each other
    inset_state.normalize_inset_area(
      cart_info.cart_initial_total_target_area());
  }

  // Shift insets so that they do not overlap
  cart_info.shift_insets_to_target_position();

  // Output to GeoJSON
  cart_info.write_geojson(
    geo_file_name,
    map_name + "_cartogram.geojson",
    output_to_stdout);

  // Stop of main function time
  time_tracker.stop("Total Time");

  // Print summary report
  time_tracker.print_summary_report();

  return EXIT_SUCCESS;
}