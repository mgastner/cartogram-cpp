#include "cartogram_info.h"
#include "constants.h"
#include "parse_arguments.h"
#include <chrono>
#include <iostream>

// Cpp Chrono for timing
typedef std::chrono::steady_clock::time_point time_point;
typedef std::chrono::steady_clock clock_time;
typedef clock_time::duration duration;
typedef std::chrono::milliseconds ms;

template <typename T> std::chrono::milliseconds inMilliseconds(T duration)
{
  return std::chrono::duration_cast<std::chrono::milliseconds>(duration);
}

int main(const int argc, const char *argv[])
{
  // Start of main function time
  const time_point start_main = clock_time::now();
  std::string geo_file_name, visual_file_name;

  // Default number of grid cells along longer Cartesian coordinate axis
  unsigned int max_n_grid_rows_or_cols;

  // Target number of points to retain after simplification
  unsigned int target_points_per_inset;
  bool world;  // World maps need special projections

  // If `triangulation` is true, we apply a cartogram projection method based
  // on the triangulation of graticule cells. It can eliminate intersections
  // that occur when the projected graticule lines are strongly curved. Only
  // use this method if the tracer points are an FTReal2d data structure.
  bool triangulation;
  bool simplify;  // Should the polygons be simplified?

  // Other boolean values that are needed to parse the command line arguments
  bool make_csv, output_equal_area, output_to_stdout, plot_density,
    plot_graticule, plot_intersections, plot_polygons, plot_quadtree,
    remove_tiny_polygons;

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
    plot_graticule,
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
    cart_info.read_geojson(geo_file_name, make_csv, &crs);
  } catch (const std::system_error &e) {
    std::cerr << "ERROR reading GeoJSON: " << e.what() << " (" << e.code()
              << ")" << std::endl;
    return EXIT_FAILURE;
  }
  std::cerr << "Coordinate reference system: " << crs << std::endl;

  // Progress measured on a scale from 0 (start) to 1 (end)
  double progress = 0.0;

  // Store total number of GeoDivs to monitor progress
  double total_geo_divs = cart_info.n_geo_divs();

  // Create std::map to store duration of each inset integrations
  std::map<std::string, ms> insets_integration_times;

  // Keep track of total time
  ms duration_initial_simplification = inMilliseconds(duration::zero()),
     duration_simplification = inMilliseconds(duration::zero()),
     duration_densification = inMilliseconds(duration::zero()),
     duration_flatten_density = inMilliseconds(duration::zero()),
     duration_fill_density = inMilliseconds(duration::zero()),
     duration_qtdt = inMilliseconds(duration::zero());

  // Project map and ensure that all holes are inside polygons
  for (auto &[inset_pos, inset_state] : *cart_info.ref_to_inset_states()) {

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

      // Simplification reduces the number of points used to represent the
      // GeoDivs in the inset, thereby reducing output file sizes and
      // run-times
      inset_state.simplify(target_points_per_inset);
    }
  }

  // Replace missing and zero target areas with positive values
  cart_info.replace_missing_and_zero_target_areas();

  // Project and exit
  if (output_equal_area) {

    // Normalize areas
    for (auto &[inset_pos, inset_state] : *cart_info.ref_to_inset_states()) {
      inset_state.normalize_inset_area(
        cart_info.cart_total_target_area(),
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

  // Iterate over insets
  for (auto &[inset_pos, inset_state] : *cart_info.ref_to_inset_states()) {

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
    inset_state.ref_to_rho_init()->allocate(lx, ly);
    inset_state.ref_to_rho_ft()->allocate(lx, ly);
    inset_state.make_fftw_plans_for_rho();
    inset_state.initialize_cum_proj();
    inset_state.set_area_errors();

    // Store initial inset area to calculate area drift
    inset_state.store_initial_area();

    // Normalize total target area to be equal to initial area
    inset_state.normalize_target_area();

    // Automatically color GeoDivs if no colors are provided
    if (inset_state.colors_empty()) {
      inset_state.auto_color();
    }
    if (plot_polygons) {

      // Write PNG and PS files if requested by command-line option
      std::string input_filename = inset_state.inset_name();
      if (plot_graticule) {
        input_filename += "_input_graticule";
      } else {
        input_filename += "_input";
      }
      std::cerr << "Writing " << input_filename << std::endl;
      inset_state.write_cairo_map(input_filename, plot_graticule);
    }

    // Remove tiny polygons below threshold
    if (remove_tiny_polygons) {
      inset_state.remove_tiny_polygons(min_polygon_area);
    }

    // We make the approximation that the progress towards generating the
    // cartogram is proportional to the number of GeoDivs that are in the
    // finished insets
    const double inset_max_frac = inset_state.n_geo_divs() / total_geo_divs;

    // Time for the initial simplification
    time_point start_initial_simplification = clock_time::now();
    if (simplify) {
      inset_state.simplify(target_points_per_inset);
    }
    time_point end_initial_simplification = clock_time::now();
    duration_initial_simplification += inMilliseconds(
      end_initial_simplification - start_initial_simplification);

    // Integration start time
    time_point start_integration = clock_time::now();

    // Start map integration
    while (inset_state.n_finished_integrations() < max_integrations &&
           (inset_state.max_area_error().value > max_permitted_area_error ||
            std::abs(inset_state.area_drift() - 1.0) > 0.01)) {
      if (qtdt_method) {
        time_point start_delaunay_t = clock_time::now();

        // Create the Delaunay triangulation
        inset_state.create_delaunay_t();
        time_point end_delaunay_t = clock_time::now();
        duration_qtdt += inMilliseconds(end_delaunay_t - start_delaunay_t);

        if (plot_quadtree) {
          const std::string quadtree_filename =
            inset_state.inset_name() + "_" +
            std::to_string(inset_state.n_finished_integrations()) +
            "_quadtree";
          std::cerr << "Writing " << quadtree_filename << std::endl;

          // Draw the resultant quadtree
          inset_state.write_quadtree(quadtree_filename);
        }
      }
      std::cerr << "Integration number "
                << inset_state.n_finished_integrations() << std::endl;

      // Calculate progress percentage. We assume that the maximum area
      // error is typically reduced to 1/5 of the previous value.
      const double ratio_actual_to_permitted_max_area_error =
        inset_state.max_area_error().value / max_permitted_area_error;
      const double n_predicted_integrations = std::max(
        (log(ratio_actual_to_permitted_max_area_error) / log(5)),
        1.0);

      // Blur density to speed up the numerics in flatten_density() below.
      // We slowly reduce the blur width so that the areas can reach their
      // target values.
      // TODO: whenever blur_width hits 0, the maximum area error will start
      //       increasing again and eventually lead to an invalid graticule
      //       cell error when projecting with triangulation. Investigate why.
      //       As a temporary fix, we set blur_width to be always positive,
      //       regardless of the number of integrations.
      double blur_width =
        std::pow(2.0, 5 - int(inset_state.n_finished_integrations()));
      // if (inset_state.n_finished_integrations() < max_integrations) {
      //   blur_width =
      //     std::pow(2.0, 5 - int(inset_state.n_finished_integrations()));
      // } else {
      //   blur_width = 0.0;
      // }
      std::cerr << "blur_width = " << blur_width << std::endl;

      // Track time needed for fill_with_density()
      time_point start_fill_density = clock_time::now();
      inset_state.fill_with_density(plot_density);
      time_point end_fill_density = clock_time::now();
      duration_fill_density +=
        inMilliseconds(end_fill_density - start_fill_density);
      if (blur_width > 0.0) {
        inset_state.blur_density(blur_width, plot_density);
      }
      if (plot_intersections) {
        inset_state.write_intersections_to_eps(intersections_resolution);
      }
      time_point start_flatten_density = clock_time::now();
      if (qtdt_method) {
        inset_state.flatten_density_with_node_vertices();
      } else {
        inset_state.flatten_density();
      }
      time_point end_flatten_density = clock_time::now();
      duration_flatten_density +=
        inMilliseconds(end_flatten_density - start_flatten_density);
      if (qtdt_method) {
        if (simplify) {
          time_point start_densify = clock_time::now();
          inset_state.densify_geo_divs_using_delaunay_t();
          time_point end_densify = clock_time::now();
          duration_densification +=
            inMilliseconds(end_densify - start_densify);
        }

        // Project using the Delaunay triangulation
        inset_state.project_with_delaunay_t();
      } else if (triangulation) {
        time_point start_densify = clock_time::now();

        // Choose diagonals that are inside graticule cells
        inset_state.fill_graticule_diagonals();

        // Densify map
        inset_state.densify_geo_divs();
        time_point end_densify = clock_time::now();
        duration_densification += inMilliseconds(end_densify - start_densify);

        // Project with triangulation
        inset_state.project_with_triangulation();
      } else {
        inset_state.project();
      }
      if (simplify) {
        time_point start_simplify = clock_time::now();
        inset_state.simplify(target_points_per_inset);
        time_point end_simplify = clock_time::now();
        duration_simplification +=
          inMilliseconds(end_simplify - start_simplify);
      }
      inset_state.increment_integration();

      // Print area drift information
      std::cerr << "Area drift: " << (inset_state.area_drift() - 1.0) * 100.0
                << "%" << std::endl;

      // Update area errors
      inset_state.set_area_errors();
      std::cerr << "max. area err: " << inset_state.max_area_error().value
                << ", GeoDiv: " << inset_state.max_area_error().geo_div
                << "\nProgress: "
                << progress + (inset_max_frac / n_predicted_integrations)
                << std::endl
                << std::endl;
    }

    // Integration end time
    time_point end_integration = clock_time::now();

    // Store integration time
    insets_integration_times[inset_pos] =
      inMilliseconds(end_integration - start_integration);
    progress += inset_max_frac;
    std::cerr << "Finished inset " << inset_pos << "\nProgress: " << progress
              << std::endl;
    if (plot_intersections) {
      inset_state.write_intersections_to_eps(intersections_resolution);
    }
    if (plot_polygons) {
      std::string output_filename = inset_state.inset_name();
      if (plot_graticule) {
        output_filename += "_output_graticule";
      } else {
        output_filename += "_output";
      }
      std::cerr << "Writing " << output_filename << std::endl;
      inset_state.write_cairo_map(output_filename, plot_graticule);
    }
    if (world) {
      std::string output_file_name =
        map_name + "_cartogram_in_smyth_projection.geojson";
      cart_info.write_geojson(
        geo_file_name,
        output_file_name,
        output_to_stdout);
      inset_state.revert_smyth_craster_projection();
    } else {

      // Rescale insets in correct proportion to each other
      inset_state.normalize_inset_area(cart_info.cart_total_target_area());
    }

    if (output_to_stdout) {
      if (qtdt_method) {
        inset_state.project_with_proj_sequence();
      } else {
        inset_state.fill_graticule_diagonals(true);
        inset_state.project_with_cum_proj();
      }
    }

    // Clean up after finishing all Fourier transforms for this inset
    inset_state.destroy_fftw_plans_for_rho();
    inset_state.ref_to_rho_init()->free();
    inset_state.ref_to_rho_ft()->free();
  }  // End of loop over insets

  // Shift insets so that they do not overlap
  cart_info.shift_insets_to_target_position();

  // Output to GeoJSON
  cart_info.write_geojson(
    geo_file_name,
    map_name + "_cartogram.geojson",
    output_to_stdout);

  // Store time when main() ended
  time_point end_main = clock_time::now();

  // Show time report
  std::cerr << std::endl;
  std::cerr << "********** Time Report **********" << std::endl;

  // Print integration times
  for (const auto &[inset_pos, inset_integration_time] :
       insets_integration_times) {
    std::cerr << "Integration Time for Inset " << inset_pos << ": "
              << inset_integration_time.count() << " ms" << std::endl;
  }
  if (qtdt_method) {
    std::cerr << "Quadtree-Delaunay T. Time: " << duration_qtdt.count()
              << " ms" << std::endl;
  }
  if (simplify) {
    std::cerr << "Initial Simplification Time: "
              << duration_initial_simplification.count() << " ms" << std::endl;
    std::cerr << "Simplification Time: " << duration_simplification.count()
              << " ms" << std::endl;
    std::cerr << "Densification Time: " << duration_densification.count()
              << " ms" << std::endl;
  }
  std::cerr << "Flatten Density Time: " << duration_flatten_density.count()
            << " ms" << std::endl;
  std::cerr << "Fill with Density Time: " << duration_fill_density.count()
            << " ms" << std::endl;
  std::cerr << "--------------------------------" << std::endl;
  std::cerr << "Total Time: " << inMilliseconds(end_main - start_main).count()
            << " ms" << std::endl;
  std::cerr << "*********************************" << std::endl;
  return EXIT_SUCCESS;
}
