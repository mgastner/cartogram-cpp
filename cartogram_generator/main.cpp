// TO DO: positional matching of argument flags

#include "constants.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include "albers_projection.h"
#include "auto_color.h"
#include "blur_density.h"
#include "check_topology.h"
#include "fill_with_density.h"
#include "flatten_density.h"
#include "project.h"
#include "read_csv.h"
#include "read_geojson.h"
#include "rescale_map.h"
#include "write_eps.h"
#include "write_to_json.h"
#include <boost/program_options.hpp>
#include <iostream>

// Functions that are called if the corresponding command-line options are
// present
void on_geometry(const std::string geometry_file_name)
{
  std::cerr << "Using geometry from file " << geometry_file_name << std::endl;
  return;
}

void on_visual_variable_file(const std::string visual_file_name)
{
  std::cerr << "Using visual variables from file "
            << visual_file_name
            << std::endl;
  return;
}

int main(const int argc, const char *argv[])
{
  using namespace boost::program_options;
  std::string geo_file_name = "", visual_file_name = ""; // Default values

  // Default number of grid cells along longer Cartesian coordinate axis.
  unsigned int long_grid_side_length = default_long_grid_side_length;

  // World maps need special projections. By default, we assume that the
  // input map is not a world map.
  bool world;

  // Other boolean values that are needed to parse the command line arguments
  bool make_csv,
       make_density_eps,
       make_polygon_eps,
       output_equal_area,
       output_to_stdout;

  // Parse command-line options. See
  // https://theboostcpplibraries.com/boost.program_options
  variables_map vm;
  try {
    options_description desc{"Options"};
    desc.add_options()(
      "help,h", "Help screen"
      )(
      "geometry,g",
      value<std::string>(&geo_file_name)
      ->required()
      ->notifier(on_geometry),
      "GeoJSON file"
      )(
      "visual_variable_file,v",
      value<std::string>(&visual_file_name)
      ->notifier(on_visual_variable_file),
      "CSV file with ID, area, and (optionally) colour"
      )(
      "output_to_stdout,s",
      value<bool>(&output_to_stdout)
      ->default_value(false)
      ->implicit_value(true),
      "Output GeoJSON to stdout"
      )(
      "output_equal_area,q",
      value<bool>(&output_equal_area)
      ->default_value(false)
      ->implicit_value(true),
      "Output equal area GeoJSON"
      )(
      "make_csv,m",
      value<bool>(&make_csv)
      ->default_value(false)
      ->implicit_value(true),
      "Boolean: create CSV file from the GeoJSON passed to the -g flag?"
      )(
      "id,i",
      value<std::string>(),
      "Column name for IDs of geographic divisions (default: 1st CSV column)"
      )(
      "area,a",
      value<std::string>(),
      "Column name for target areas (default: 2nd CSV column)"
      )(
      "color,c",
      value<std::string>(),
      "Column name for colors (default: \"Color\" or \"Colour\")"
      )(
      "inset,n",
      value<std::string>(),
      "Column name for insets (default: \"Inset\")"
      )(
      "long_grid_side_length,l",
      value<unsigned int>(&long_grid_side_length),
      "Number of grid cells along longer Cartesian coordinate axis"
      )(
      "world,w",
      value<bool>(&world)
      ->default_value(false)
      ->implicit_value(true),
      "Boolean: is input a world map in longitude-latitude format?"
      )(
      "polygons_to_eps,e",
      value<bool>(&make_polygon_eps)
      ->default_value(false)
      ->implicit_value(true),
      "Boolean: make EPS image of input and output?"
      )(
      "density_to_eps,d",
      value<bool>(&make_density_eps)
      ->default_value(false)
      ->implicit_value(true),
      "Boolean: make EPS images *_density_*.eps?"
      );
    store(parse_command_line(argc, argv, desc), vm);
    if (vm.count("help") || argc == 1) {
      std::cerr << desc << '\n';
      return EXIT_SUCCESS;
    } else {
      notify(vm);  // Triggers notifier functions such as on_geometry()
    }
  } catch (const error &ex) {
    std::cerr << "ERROR: " << ex.what() << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize cart_info. It contains all information about the cartogram
  // that needs to be handled by functions called from main().
  CartogramInfo cart_info(world, visual_file_name, make_density_eps);
  if (!make_csv) {

    // Read visual variables (e.g. area, color) from CSV
    try {
      read_csv(vm, &cart_info);
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

  // Read geometry
  try {
    read_geojson(geo_file_name, &cart_info, make_csv);
  } catch (const std::system_error& e) {
    std::cerr << "ERROR: "
              << e.what()
              << " ("
              << e.code()
              << ")"
              << std::endl;
    return EXIT_FAILURE;
  }

  // Determine name of input map
  std::string map_name = geo_file_name;
  if (map_name.find_last_of("/\\") != std::string::npos) {
    map_name = map_name.substr(map_name.find_last_of("/\\") + 1);
  }
  if (map_name.find('.') != std::string::npos) {
    map_name = map_name.substr(0, map_name.find('.'));
  }

  // Loop over insets
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
    CGAL::Bbox_2 bb = inset_state.bbox();
    if (bb.xmin() >= -180.0 && bb.xmax() <= 180.0 &&
        bb.ymin() >= -90.0 && bb.ymax() <= 90.0) {

      // If yes, transform the coordinates with the Albers projection
      try {
        transform_to_albers_projection(&inset_state);
      } catch (const std::system_error& e) {
        std::cerr << "ERROR: "
                  << e.what()
                  << " ("
                  << e.code()
                  << ")"
                  << std::endl;
        return EXIT_FAILURE;
      }
    } else if (output_equal_area) {
      std::cerr << "ERROR: Input GeoJSON is not a longitude-latitude map."
                << std::endl;
      return EXIT_FAILURE;
    }

    // Determine the name of the inset
    std::string inset_name = map_name;
    if (cart_info.n_insets() > 1) {
      inset_name = inset_name + "_" + inset_state.pos();
      std::cerr << "\nWorking on inset at position: "
                << inset_state.pos()
                << std::endl;
    }
    inset_state.set_inset_name(inset_name);
    if (output_equal_area) {
      normalize_inset_area(&inset_state,
                           cart_info.total_cart_target_area(),
                           output_equal_area);
    } else {

      // Rescale map to fit into a rectangular box [0, lx] * [0, ly].
      rescale_map(long_grid_side_length,
                  &inset_state,
                  cart_info.is_world_map());

      // Set up Fourier transforms
      unsigned int lx = inset_state.lx();
      unsigned int ly = inset_state.ly();
      inset_state.ref_to_rho_init()->allocate(lx, ly);
      inset_state.ref_to_rho_ft()->allocate(lx, ly);
      inset_state.make_fftw_plans_for_rho();

      // Set initial area errors
      inset_state.set_area_errors();

      // Fill density to fill horizontal adjacency map
      fill_with_density(&inset_state,
                        cart_info.trigger_write_density_to_eps());

      // Automatically color GeoDivs if no colors are provided
      if (inset_state.colors_empty()) {
        auto_color(&inset_state);
      }

      // Write EPS if requested by command-line option
      if (make_polygon_eps) {
        std::cerr << "Writing " << inset_name << "_input.eps" << std::endl;
        write_map_to_eps((inset_name + "_input.eps"), &inset_state);
      }

      // Start map integration
      while (inset_state.n_finished_integrations() < max_integrations &&
             inset_state.max_area_error() > max_permitted_area_error) {
        std::cerr << "Integration number "
                  << inset_state.n_finished_integrations()
                  << std::endl;

        // TODO: THIS IF-CONDITION IS INELEGANT
        if (inset_state.n_finished_integrations()  >  0) {
          fill_with_density(&inset_state,
                            cart_info.trigger_write_density_to_eps());
        }
        if (inset_state.n_finished_integrations() == 0) {
          blur_density(5.0,
                       &inset_state,
                       cart_info.trigger_write_density_to_eps());
        } else {
          blur_density(0.0,
                       &inset_state,
                       cart_info.trigger_write_density_to_eps());
        }
        flatten_density(&inset_state);
        project(&inset_state);
        inset_state.increment_integration();

        // Update area errors
        inset_state.set_area_errors();
      }

      // Print EPS of cartogram
      if (make_polygon_eps) {
        std::cerr << "Writing "
                  << inset_state.inset_name()
                  << "_output.eps" << std::endl;
        write_map_to_eps((inset_state.inset_name() + "_output.eps"),
                         &inset_state);
      }

      // Rescale insets in correct proportion to each other
      normalize_inset_area(&inset_state,
                           cart_info.total_cart_target_area());

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
  nlohmann::json cart_json = cgal_to_json(&cart_info);
  write_to_json(cart_json,
                geo_file_name,
                output_file_name,
                std::cout,
                output_to_stdout);
  return EXIT_SUCCESS;
}
