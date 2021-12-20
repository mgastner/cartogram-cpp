#include "albers_projection.h"
#include "auto_color.h"
#include "blur_density.h"
#include "cartogram_info.h"
#include "check_topology.h"
#include "constants.h"
#include "densification_points.h"
#include "densify.h"
#include "fill_with_density.h"
#include "flatten_density.h"
#include "geo_div.h"
#include "inset_state.h"
#include "project.h"
#include "read_csv.h"
#include "read_geojson.h"
#include "rescale_map.h"
#include "simplify_inset.h"
#include "write_eps.h"
#include "write_geojson.h"
#include "xy_point.h"
#include "parse_arguments.h"
#include <iostream>
#include <cmath>

int main(const int argc, const char *argv[])
{
  std::string geo_file_name, visual_file_name; // Default values

  // Default number of grid cells along longer Cartesian coordinate axis.
  unsigned int long_grid_side_length = default_long_grid_side_length;

  // World maps need special projections. By default, we assume that the
  // input map is not a world map.
  bool world;

  // Another cartogram projection method based on triangulation of graticule
  // cells. It can reduce or even eliminate intersections that occur when
  // projecting "naively".
  bool triangulation;

  // Shall the polygons be simplified?
  bool simplify;

  // Other boolean values that are needed to parse the command line arguments
  bool make_csv,
       make_polygon_eps,
       output_equal_area,
       output_to_stdout,
       plot_density,
       plot_graticule;

  // Parse command-line arguments
  argparse::ArgumentParser arguments = parsed_arguments(
    argc,
    argv,
    geo_file_name,
    visual_file_name,
    long_grid_side_length,
    world,
    triangulation,
    simplify,
    make_csv,
    make_polygon_eps,
    output_equal_area,
    output_to_stdout,
    plot_density,
    plot_graticule);


  // Initialize cart_info. It contains all information about the cartogram
  // that needs to be handled by functions called from main().
  CartogramInfo cart_info(world, visual_file_name);
  if (!make_csv) {

    // Read visual variables (e.g. area, color) from CSV
    try {
      read_csv(arguments, &cart_info);
    } catch (const std::system_error& e) {
      std::cerr << "ERROR: "
                << e.what()
                << " ("
                << e.code()
                << ")"
                << std::endl;
      return EXIT_FAILURE;
    } catch (const std::runtime_error& e) {

      // Likely due to invalid CSV file
      std::cerr << "ERROR: "
                << e.what()
                << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Read geometry. If the GeoJSON does not explicitly contain a "crs" field,
  // assume that the coordinates are in longitude and latitude.
  std::string crs = "+proj=longlat";
  try {
    read_geojson(geo_file_name, make_csv, &crs, &cart_info);
  } catch (const std::system_error& e) {
    std::cerr << "ERROR: "
              << e.what()
              << " ("
              << e.code()
              << ")"
              << std::endl;
    return EXIT_FAILURE;
  }
  std::cerr << "Coordinate reference system: " << crs << std::endl;

  // Progress percentage
  double progress = 0.0;

  // Store total number of GeoDivs to monitor progress
  double total_geo_divs = 0;
  for (const auto &inset_info : *cart_info.ref_to_inset_states()) {

    // 'auto' will automatically deduce the const qualifier.
    auto &inset_state = inset_info.second;
    total_geo_divs += inset_state.n_geo_divs();
  }

  // Determine name of input map
  std::string map_name = geo_file_name;
  if (map_name.find_last_of("/\\") != std::string::npos) {
    map_name = map_name.substr(map_name.find_last_of("/\\") + 1);
  }
  if (map_name.find('.') != std::string::npos) {
    map_name = map_name.substr(0, map_name.find('.'));
  }

  // Project map and ensure that all holes are inside polygons
  for (auto &[inset_pos, inset_state] : *cart_info.ref_to_inset_states()) {

    // Check for errors in the input topology
    try {
      holes_inside_polygons(&inset_state);
    } catch (const std::system_error& e) {
      std::cerr << "ERROR: "
                << e.what()
                << " ("
                << e.code()
                << ")"
                << std::endl;
      return EXIT_FAILURE;
    }

    // Can the coordinates be interpreted as longitude and latitude?
    const Bbox bb = inset_state.bbox();
    if (bb.xmin() >= -180.0 && bb.xmax() <= 180.0 &&
        bb.ymin() >= -90.0 && bb.ymax() <= 90.0 &&
        crs == "+proj=longlat") {

      // If yes, transform the coordinates with the Albers projection
      transform_to_albers_projection(&inset_state);
    } else if (output_equal_area) {
      std::cerr << "ERROR: Input GeoJSON is not a longitude-latitude map."
                << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Replace missing and zero target areas with absolute values
  cart_info.replace_missing_and_zero_target_areas();

  // Loop over insets
  for (auto &[inset_pos, inset_state] : *cart_info.ref_to_inset_states()) {

    // Determine the name of the inset
    std::string inset_name = map_name;
    if (cart_info.n_insets() > 1) {
      inset_name = inset_name + "_" + inset_pos;
      std::cerr << "\nWorking on inset at position: "
                << inset_pos
                << std::endl;
    }
    inset_state.set_inset_name(inset_name);

    // Simplify inset if -s flag is passed
    if (simplify) {
      simplify_inset(&inset_state);
    }
    if (output_equal_area) {
      normalize_inset_area(&inset_state,
                           cart_info.cart_total_target_area(),
                           output_equal_area);
    } else {

      // Rescale map to fit into a rectangular box [0, lx] * [0, ly].
      rescale_map(long_grid_side_length,
                  &inset_state,
                  cart_info.is_world_map());

      // Set up Fourier transforms
      const unsigned int lx = inset_state.lx();
      const unsigned int ly = inset_state.ly();
      inset_state.ref_to_rho_init()->allocate(lx, ly);
      inset_state.ref_to_rho_ft()->allocate(lx, ly);
      inset_state.make_fftw_plans_for_rho();
      inset_state.initialize_cum_proj();

      // Set initial area errors
      inset_state.set_area_errors();

      // TODO: AN UPDATED VERSION OF auto_color() SHOULD AVOID HAVING TO
      // CALL fill_with_density() TO GET ADJACENCY.

      // Fill density to fill horizontal adjacency map
      fill_with_density(plot_density, &inset_state);

      // Automatically color GeoDivs if no colors are provided
      if (inset_state.colors_empty()) {
        auto_color(&inset_state);
      }

      // Write EPS if requested by command-line option
      if (make_polygon_eps) {
        std::string eps_input_filename = inset_state.inset_name();
        if (plot_graticule) {
          eps_input_filename += "_input_graticule.eps";
        } else {
          eps_input_filename += "_input.eps";
        }
        std::cerr << "Writing " << eps_input_filename << std::endl;
        write_map_to_eps(eps_input_filename, plot_graticule, &inset_state);
      }

      // We make the approximation that the progress towards generating the
      // cartogram is proportional to the number of GeoDivs that are in the
      // finished insets
      const double inset_max_frac = inset_state.n_geo_divs() / total_geo_divs;


      // Start map integration
      while (inset_state.n_finished_integrations() < max_integrations &&
             inset_state.max_area_error().value > max_permitted_area_error) {
        std::cerr << "Integration number "
                  << inset_state.n_finished_integrations()
                  << std::endl;

        // Calculate progress percentage. We assume that the maximum area
        // error is typically reduced to 1/5 of the previous value.
        const double ratio_actual_to_permitted_max_area_error =
          inset_state.max_area_error().value / max_permitted_area_error;
        const double n_predicted_integrations =
          std::max((log(ratio_actual_to_permitted_max_area_error) / log(5)),
                   1.0);
        double blur_width;
        if (inset_state.n_finished_integrations() < 10) {
          blur_width =
            std::pow(2.0, 5 - int(inset_state.n_finished_integrations()));
        } else {
          blur_width = 0.0;
        }
        std::cerr << "blur_width = " << blur_width << std::endl;

        // TODO: THIS if-CONDITION IS INELEGANT. IN AN UPDATED VERSION OF
        // auto_color() THE EARLIER fill_with_density() SHOULD BE REMOVED
        // AND THE if-CONDITION BE DROPPED.
        if (inset_state.n_finished_integrations() > 0) {
          fill_with_density(plot_density, &inset_state);
        }

        // TODO: DO NOT CARRY OUT blur_density() IF blur_width=0.0
        blur_density(blur_width, plot_density, &inset_state);
        flatten_density(&inset_state);
        if (triangulation) {

          // Choose diagonals that are inside graticule cells
          fill_graticule_diagonals(&inset_state);

          // Densify map
          // TODO: It would make sense to turn densified_geo_divs() into a
          // method of InsetState. Then the next command could be written more
          // simply as inset_state.densify_geo_divs().
          inset_state.set_geo_divs(
            densified_geo_divs(inset_state.geo_divs()));

          // Project with triangulation
          project_with_triangulation(&inset_state);
        } else {
          project(&inset_state);
        }
        if (simplify) {
          simplify_inset(&inset_state);
        }
        inset_state.increment_integration();

        // Update area errors
        inset_state.set_area_errors();
        std::cerr << "max. area err: "
                  << inset_state.max_area_error().value
                  << ", GeoDiv: "
                  << inset_state.max_area_error().geo_div
                  << std::endl;
        std::cerr << "Progress: "
                  << progress + (inset_max_frac / n_predicted_integrations)
                  << std::endl
                  << std::endl;
      }
      progress += inset_max_frac;
      std::cerr << "Finished inset "
                << inset_pos
                << "\nProgress: "
                << progress
                << std::endl;

      // Print EPS of cartogram
      if (make_polygon_eps) {
        std::string eps_output_filename = inset_state.inset_name();
        if (plot_graticule) {
          eps_output_filename += "_output_graticule.eps";
        } else {
          eps_output_filename += "_output.eps";
        }
        std::cerr << "Writing "
                  << eps_output_filename << std::endl;
        write_map_to_eps(eps_output_filename, plot_graticule,
                         &inset_state);
      }

      // Rescale insets in correct proportion to each other
      normalize_inset_area(&inset_state,
                           cart_info.cart_total_target_area());

      // Clean up after finishing all Fourier transforms for this inset
      inset_state.destroy_fftw_plans_for_rho();
      inset_state.ref_to_rho_init()->free();
      inset_state.ref_to_rho_ft()->free();
    } // End of loop over insets
  }

  // Shift insets so that they do not overlap
  shift_insets_to_target_position(&cart_info);

  // Output to GeoJSON
  std::string output_file_name;
  if (output_equal_area) {
    output_file_name = map_name + "_equal_area.geojson";
  } else {
    output_file_name = map_name + "_cartogram.geojson";
  }
  const nlohmann::json cart_json = cgal_to_json(&cart_info);
  write_geojson(cart_json,
                geo_file_name,
                output_file_name,
                std::cerr,
                output_to_stdout,
                &cart_info);
  return EXIT_SUCCESS;
}
