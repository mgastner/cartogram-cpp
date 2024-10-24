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
  unsigned int long_grid_side_length;

  // Target number of points to retain after simplification
  unsigned int target_points_per_inset;
  bool world;  // World maps need special projections

  // If `triangulation` is true, we apply a cartogram projection method based
  // on the triangulation of grid cells. It can eliminate intersections
  // that occur when the projected grid lines are strongly curved. Only
  // use this method if the tracer points are an FTReal2d data structure.
  bool triangulation;
  bool simplify;  // Should the polygons be simplified?

  // If `rays` is true, we use the ray-shooting method to fill the grid cells.
  bool rays;

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
    long_grid_side_length,
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
    plot_quadtree,
    rays);

  // Initialize cart_info. It contains all the information about the cartogram
  // that needs to be handled by functions called from main().
  CartogramInfo cart_info(world, visual_file_name);

  // Determine name of input map based on the geo_file_name and store it
  std::string map_name = cart_info.set_map_name(geo_file_name);
  if (!make_csv) {

    // Read visual variables (e.g., area and color) from CSV
    cart_info.read_csv(arguments);
  }

  // Read geometry. If the GeoJSON does not explicitly contain a "crs" field,
  // we assume that the coordinates are in longitude and latitude.
  std::string crs = "+proj=longlat";
  cart_info.read_geojson(geo_file_name, make_csv, crs);
  std::cerr << "Coordinate reference system: " << crs << std::endl;

  // Store total number of GeoDivs to monitor progress
  double total_geo_divs = cart_info.n_geo_divs();

  // Project map and ensure that all holes are inside polygons
  for (auto &[inset_pos, inset_state] : cart_info.ref_to_inset_states()) {

    // Start of inset time
    time_tracker.start("Inset " + inset_pos);

    // Check for errors in the input topology
    inset_state.check_topology();

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
      map_name + "_equal_area",
      false);
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
    inset_state.rescale_map(long_grid_side_length, cart_info.is_world_map());

    if (output_to_stdout) {

      // Store original coordinates
      inset_state.store_original_geo_divs();
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

    // Output rescaled GeoJSON
    cart_info.write_geojson(
      geo_file_name,
      // processed = simplified + rescaled
      // and potentially projected + small polygons removed
      map_name + "_input_processed",
      false);

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

    // Store initial inset area to calculate area drift,
    // set area errors based on this initial_area
    inset_state.store_initial_area();
    inset_state.set_area_errors();

    // Store initial target area to normalize inset areas
    inset_state.store_initial_target_area();

    // Normalize total target area to be equal to initial area
    inset_state.normalize_target_area();

    // Automatically color GeoDivs if no colors are provided
    if (inset_state.colors_empty()) {
      inset_state.auto_color();
    }
    if (plot_polygons) {

      // Write input of SVG files if requested by command-line option
      inset_state.write_cairo_map(
        inset_state.inset_name() + "_input",
        plot_grid);
    }

    // Remove tiny polygons below threshold
    if (remove_tiny_polygons) {
      inset_state.remove_tiny_polygons(min_polygon_area);
    }

    time_tracker.start("Integration Inset " + inset_pos);

    // Print initial values, and initial progress (0)
    progress_tracker.print_progress_mid_integration(inset_state);

    // Start map integration
    while (inset_state.continue_integrating()) {
      // File prefix for output files for this integration
      const std::string file_prefix =
        inset_state.inset_name() + "_" +
        std::to_string(inset_state.n_finished_integrations());

      if (qtdt_method) {

        // Create the Delaunay triangulation
        time_tracker.start("Delaunay Triangulation");
        inset_state.create_and_store_quadtree_cell_corners();
        time_tracker.stop("Delaunay Triangulation");

        time_tracker.start("Delaunay Triangulation");
        inset_state.create_delaunay_t();
        time_tracker.stop("Delaunay Triangulation");

        if (plot_quadtree) {

          // Draw the resultant quadtree and Delaunay triangulation
          inset_state.write_quadtree(file_prefix + "_quadtree");
          inset_state.write_delaunay_triangles(
            file_prefix + "_delaunay_t",
            false);
        }
      }

      if (rays) {
        // Fill density using ray-shooting method
        time_tracker.start("Fill with Density (Ray Shooting Method)");
        inset_state.fill_with_density_rays(plot_density);
        time_tracker.stop("Fill with Density (Ray Shooting Method)");
      } else {
        time_tracker.start("Fill with Density (Clipping Method)");
        inset_state.fill_with_density_clip(plot_density);
        time_tracker.stop("Fill with Density (Clipping Method)");
      }

      const double blur_width = inset_state.blur_width();
      if (blur_width > 0.0) {
        time_tracker.start("Blur");
        inset_state.blur_density(blur_width, plot_density);
        time_tracker.stop("Blur");
      }
      if (qtdt_method) {

        time_tracker.start("Flatten Density (Quadtree Method)");
        if (!inset_state.flatten_density_with_node_vertices()) {

          // Flatten density has failed. Incrrease blur width and try again
          time_tracker.stop("Flatten Density (Quadtree Method)");
          inset_state.increment_n_fails_during_flatten_density();
          continue;
        }

        // Flatten density passed
        time_tracker.stop("Flatten Density (Quadtree Method)");

        if (plot_quadtree) {
          inset_state.write_delaunay_triangles(
            file_prefix + "_delaunay_t_after_flatten",
            true);
        }
      } else {

        // Using entire grid
        time_tracker.start("Flatten Density (Full Grid Method)");
        inset_state.flatten_density();
        time_tracker.stop("Flatten Density (Full Grid Method)");
      }

      if (qtdt_method) {
        time_tracker.start("Update Delanuay Triangulation");
        inset_state.update_delaunay_t();
        time_tracker.stop("Update Delanuay Triangulation");

        if (plot_quadtree) {
          inset_state.write_delaunay_triangles(
            file_prefix + "_updated_delaunay_t_projected",
            true);
        }

        if (simplify) {

          time_tracker.start("Densification (using Delanuay Triangles)");
          inset_state.densify_geo_divs_using_delaunay_t();
          time_tracker.stop("Densification (using Delanuay Triangles)");
        }

        // Project using the Delaunay triangulation
        time_tracker.start("Project (Delanuay Triangulation)");
        inset_state.project_with_delaunay_t(output_to_stdout);
        time_tracker.stop("Project (Delanuay Triangulation)");
      } else if (triangulation) {

        // Choose diagonals that are inside grid cells, then densify.
        time_tracker.start("Densification (using Grid Diagonals)");
        inset_state.fill_grid_diagonals();
        inset_state.densify_geo_divs();
        time_tracker.stop("Densification (using Grid Diagonals)");

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
      if (simplify) {

        time_tracker.start("Simplification");
        inset_state.simplify(target_points_per_inset);
        time_tracker.stop("Simplification");
      }
      if (plot_intersections) {
        inset_state.write_intersections_image(intersections_resolution);
      }

      // Update area errors
      inset_state.set_area_errors();
      inset_state.adjust_grid();
      progress_tracker.print_progress_mid_integration(inset_state);
      inset_state.increment_integration();
    }
    time_tracker.stop("Integration Inset " + inset_pos);

    // Update and display progress information
    std::cerr << "Finished inset " << inset_pos << std::endl;
    progress_tracker.update_and_print_progress_end_integration(inset_state);

    if (plot_polygons) {
      inset_state.write_cairo_map(
        inset_state.inset_name() + "_output",
        plot_grid);
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
    map_name + "_cartogram",
    output_to_stdout);

  // Stop of main function time
  time_tracker.stop("Total Time");

  // Print summary report
  time_tracker.print_summary_report();

  return EXIT_SUCCESS;
}