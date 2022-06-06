#include "cartogram_info.h"
#include "constants.h"
#include "parse_arguments.h"

#include <chrono>
#include <cmath>
#include <iostream>

typedef std::chrono::steady_clock::time_point time_point;
typedef std::chrono::steady_clock clock_time;
typedef std::chrono::milliseconds ms;

template <typename T> std::chrono::milliseconds inMilliseconds(T duration)
{
  return std::chrono::duration_cast<std::chrono::milliseconds>(duration);
}

int main(const int argc, const char *argv[])
{

  // Start of main function time
  time_point start_main = clock_time::now();
  std::string geo_file_name, visual_file_name;

  // Default number of grid cells along longer Cartesian coordinate axis
  unsigned int long_grid_length = default_long_grid_length;

  // Target number of points to retain after simplification
  unsigned int target_points_per_inset = default_target_points_per_inset;
  bool world;  // World maps need special projections

  // Another cartogram projection method based on triangulation of grid
  // cells. It can eliminate intersections that occur when the projected
  // grid lines are strongly curved.
  bool triangulation;

  // Shall the polygons be simplified?
  bool simplify;

  // Other boolean values that are needed to parse the command line arguments

  bool make_csv, produce_map_image, image_format_ps, output_equal_area,
    output_to_stdout, plot_density, plot_grid, plot_grid_heatmap,
    plot_intersections, crop;

  // Parse command-line arguments
  argparse::ArgumentParser arguments = parsed_arguments(
    argc,
    argv,
    geo_file_name,
    visual_file_name,
    long_grid_length,
    target_points_per_inset,
    world,
    triangulation,
    simplify,
    make_csv,
    produce_map_image,
    image_format_ps,
    output_equal_area,
    output_to_stdout,
    plot_density,
    plot_grid,
    plot_grid_heatmap,
    plot_intersections,
    crop);

  // Initialize cart_info. It contains all information about the cartogram
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
  std::cout << geo_file_name << "\n";

  cart_info.set_map_name(map_name);

  // Parsing start time
  time_point start_parse = clock_time::now();

  if (!make_csv) {
    // Read visual variables (e.g. area, color) from CSV
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

  // Parsing end time
  time_point end_parse = clock_time::now();

  std::cerr << "Coordinate reference system: " << crs << std::endl;

  // Progress measured on a scale from 0 (start) to 1 (end)
  double progress = 0.0;

  // Store total number of GeoDivs to monitor progress
  double total_geo_divs = cart_info.n_geo_divs();

  // Projection and simplification start time
  time_point start_proj_sim = clock_time::now();

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
    // TODO: the "crs" field for GeoJSON files seems to be deprecated. However,
    // in earlier specifications, the coordinate reference system used to be
    // written in the format specified here:
    // https://geojson.org/geojson-spec.html#coordinate-reference-system-objects.
    // It may be a good idea to make a list of possible entries corresponding
    // to longitude and lattitude projection. "urn:ogc:def:crs:OGC:1.3:CRS84"
    // is one such entry.
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
    } else if (output_equal_area || plot_grid_heatmap) {
      std::cerr << "ERROR: Input GeoJSON is not a longitude-latitude map."
                << std::endl;
      return EXIT_FAILURE;
    }
    if (simplify) {

      // Simplification reduces the number of points used to represent the
      // GeoDivs in the inset, thereby reducing output file sizes and run
      // times
      inset_state.simplify(target_points_per_inset);
    }
  }

  // Albers projection end time
  time_point end_proj_sim = clock_time::now();

  // Replace missing and zero target areas with absolute values
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
      std::cout,
      output_to_stdout);
    return EXIT_SUCCESS;
  }

  // Create map to store duration of each inset integrations
  std::map<std::string, ms> insets_integration_times;

  // Iterate over insets
  for (auto &[inset_pos, inset_state] : *cart_info.ref_to_inset_states()) {
    // Determine the name of the inset
    std::string inset_name = map_name;
    if (cart_info.n_insets() > 1) {
      inset_name = inset_name + "_" + inset_pos;
      std::cerr << "\nWorking on inset at position: " << inset_pos
                << std::endl;
    }
    inset_state.set_inset_name(inset_name);

    // Rescale map to fit into a rectangular box [0, lx] * [0, ly]
    inset_state.rescale_map(long_grid_length, cart_info.is_world_map());

    // Store original coordinates
    inset_state.store_original_geo_divs();

    // Set up Fourier transforms
    const unsigned int lx = inset_state.lx();
    const unsigned int ly = inset_state.ly();
    inset_state.ref_to_rho_init()->allocate(lx, ly);
    inset_state.ref_to_rho_ft()->allocate(lx, ly);
    inset_state.make_fftw_plans_for_rho();
    inset_state.initialize_cum_proj();
    inset_state.initialize_identity_proj();
    inset_state.set_area_errors();

    // Automatically color GeoDivs if no colors are provided
    if (inset_state.colors_empty()) {
      inset_state.auto_color();
    }

    // Write PNG and PS files if requested by command-line option
    if (produce_map_image) {
      std::string input_filename = inset_state.inset_name();

      // Update filename with grid info
      plot_grid ? input_filename += "_input_grid" : input_filename += "_input";

      // Update extension
      image_format_ps ? input_filename += ".ps" : input_filename += ".svg";

      std::cerr << "Writing " << input_filename << std::endl;
      inset_state
        .write_map_image(input_filename, true, plot_grid, image_format_ps);
    }

    // We make the approximation that the progress towards generating the
    // cartogram is proportional to the number of GeoDivs that are in the
    // finished insets
    const double inset_max_frac = inset_state.n_geo_divs() / total_geo_divs;

    // Integration start time
    time_point start_integration = clock_time::now();

    // Start map integration
    while (inset_state.n_finished_integrations() < max_integrations &&
           inset_state.max_area_error().value > max_permitted_area_error) {
      std::cerr << "Integration number "
                << inset_state.n_finished_integrations() << std::endl;

      // Blur density to speed up the numerics in flatten_density() below.
      // We slowly reduce the blur width so that the areas can reach their
      // target values.
      // TODO: whenever blur_width hits 0, the maximum area error will start
      // increasing again and eventually lead to an invalid grid cell
      // error when projecting with triangulation. Investigate why. As a
      // temporary fix, we set blur_width to be always positive, regardless
      // of the number of integrations. This error also seems to occur when
      // the starting blur width is too low to begin with.
      // TODO: Add option to customise starting blur width.
      double blur_width =
        std::pow(2.0, 5 - int(inset_state.n_finished_integrations()));
      // if (inset_state.n_finished_integrations() < max_integrations) {
      //   blur_width =
      //     std::pow(2.0, 3 - int(inset_state.n_finished_integrations()));
      // } else {
      //   blur_width = 0.0;
      // }
      std::cerr << "blur_width = " << blur_width << std::endl;

      inset_state.fill_with_density(
        plot_density,
        plot_grid_heatmap,
        image_format_ps);
      if (blur_width > 0.0) {
        inset_state.blur_density(blur_width, plot_density, image_format_ps);
      }
      if (plot_intersections) {
        inset_state.write_intersections_image(
          intersections_resolution,
          image_format_ps);
      }
      inset_state.flatten_density();
      if (triangulation) {

        // Choose diagonals that are inside grid cells
        inset_state.fill_grid_diagonals();

        // Densify map
        inset_state.densify_geo_divs();

        // Project with triangulation
        inset_state.project_with_triangulation();
      } else {
        inset_state.project();
      }
      if (simplify) {
        inset_state.simplify(target_points_per_inset);
      }
      inset_state.increment_integration();

      // Update area errors
      inset_state.set_area_errors();

      // Calculate progress percentage. We assume that the maximum area
      // error is typically reduced to 1/10 of the previous value.
      const double ratio_actual_to_permitted_max_area_error = std::max(
        inset_state.max_area_error().value / max_permitted_area_error,
        1.0);
      const double n_predicted_integrations =
        std::max((log10(ratio_actual_to_permitted_max_area_error)), 1.0);

      std::cerr << "max. area err: " << inset_state.max_area_error().value
                << ", GeoDiv: " << inset_state.max_area_error().geo_div
                << "\nProgress: "
                << progress + (inset_max_frac / n_predicted_integrations)
                << std::endl
                << std::endl;
    }

    // Integration end time
    time_point end_integration = clock_time::now();

    // Add integration time to the map
    insets_integration_times[inset_pos] =
      inMilliseconds(end_integration - start_integration);
    progress += inset_max_frac;
    std::cerr << "Finished inset " << inset_pos << "\nProgress: " << progress
              << std::endl;

    if (plot_intersections) {
      inset_state.write_intersections_image(
        intersections_resolution,
        image_format_ps);
    }

    // Print PS files of cartogram
    if (produce_map_image) {
      std::string output_filename = inset_state.inset_name();

      // Update filename with grid info
      plot_grid ? output_filename += "_output_grid"
                : output_filename += "_output";

      // Update extension
      image_format_ps ? output_filename += ".ps" : output_filename += ".svg";

      std::cerr << "Writing " << output_filename << std::endl;
      inset_state
        .write_map_image(output_filename, true, plot_grid, image_format_ps);
    }

    if (plot_grid_heatmap) {
      // Produce equal-area grid heatmap
      std::string output_filename =
        inset_state.inset_name() + "_equal_area_grid_heatmap";

      // Update extension
      image_format_ps ? output_filename += ".ps" : output_filename += ".svg";
      std::cerr << "Writing " << output_filename << std::endl;
      inset_state.write_grid_heatmap_image(
        output_filename,
        true,
        image_format_ps,
        crop);

      // Produce cartogram grid heatmap
      output_filename = inset_state.inset_name() + "_cartogram_grid_heatmap";

      // Update extension
      image_format_ps ? output_filename += ".ps" : output_filename += ".svg";

      std::cerr << "Writing " << output_filename << std::endl;
      inset_state.write_grid_heatmap_image(
        output_filename,
        false,
        image_format_ps,
        crop);
    }
    if (world) {
      std::string output_file_name =
        map_name + "_cartogram_in_smyth_projection.geojson";
      cart_info.write_geojson(
        geo_file_name,
        output_file_name,
        std::cout,
        output_to_stdout);
      inset_state.revert_smyth_craster_projection();
    } else {
      // Rescale insets in correct proportion to each other
      inset_state.normalize_inset_area(cart_info.cart_total_target_area());
    }

    // Clean up after finishing all Fourier transforms for this inset
    inset_state.destroy_fftw_plans_for_rho();
    inset_state.ref_to_rho_init()->free();
    inset_state.ref_to_rho_ft()->free();
  }
  // End of loop over insets

  // Output a density heatmap's bar
  if (plot_density) {
    std::string output_filename = "density_bar";

    // Update extension
    image_format_ps ? output_filename += ".ps" : output_filename += ".svg";
    std::cerr << "Writing " << output_filename << std::endl;
    // write_density_bar_image(output_filename, image_format_ps);
  }

  // Shift insets so that they do not overlap
  cart_info.shift_insets_to_target_position();

  // Output to GeoJSON
  cart_info.write_geojson(
    geo_file_name,
    map_name + "_cartogram.geojson",
    std::cout,
    output_to_stdout);

  // End of main function time
  time_point end_main = clock_time::now();

  // Calculate differences in time
  ms parsing_time = inMilliseconds(end_parse - start_parse);
  ms proj_sim_time = inMilliseconds(end_proj_sim - start_proj_sim);
  ms total_time = inMilliseconds(end_main - start_main);

  // Show Time Report
  std::cerr << std::endl;
  std::cerr << "********** Time Report **********" << std::endl;
  std::cerr << "Parsing time: " << parsing_time.count() << " ms" << std::endl;
  std::cerr << "Projection & Simplification (if any) time: "
            << proj_sim_time.count() << " ms" << std::endl;

  // Iterate over the map and print integration times
  for (auto [inset_pos, inset_integration_time] : insets_integration_times) {
    std::cerr << "Integration time for inset " << inset_pos << ": "
              << inset_integration_time.count() << " ms" << std::endl;
  }
  std::cerr << "Total time: " << total_time.count() << " ms" << std::endl;
  std::cerr << "*********************************" << std::endl;

  return EXIT_SUCCESS;
}
