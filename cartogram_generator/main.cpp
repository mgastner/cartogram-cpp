// TO DO: positional matching of argument flags

#include "constants.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include "blur_density.h"
#include "fill_with_density.h"
#include "flatten_density.h"
#include "project.h"
#include "read_csv.h"
#include "read_geojson.h"
#include "rescale_map.h"
#include "write_eps.h"
#include "check_topology.h"
#include "write_to_json.h"
#include "auto_color.h"
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
  int long_grid_side_length = default_long_grid_side_length;

  // World maps need special projections. By default, we assume that the
  // input map is not a world map.
  bool world;

  // Other boolean values that are needed to parse the command line arguments
  bool polygons_to_eps, density_to_eps, make_csv;

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
      "make_csv,m",
      value<bool>(&make_csv)
      ->default_value(false)
      ->implicit_value(true),
      "Boolean: create a CSV file from the GeoJSON file passed to the -g flag?"
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
      "Column name for colors (assumed column name: \"Color\" or \"Colour\")"
      )(
      "inset,i",
      value<std::string>(),
      "Column name for insets (assumed column name: \"Inset\")"
      )(
      "long_grid_side_length,l",
      value<int>(&long_grid_side_length),
      "Number of grid cells along longer Cartesian coordinate axis"
      )(
      "world,w",
      value<bool>(&world)
      ->default_value(false)
      ->implicit_value(true),
      "Boolean: is input a world map in longitude-latitude format?"
      )(
      "polygons_to_eps,e",
      value<bool>(&polygons_to_eps)
      ->default_value(false)
      ->implicit_value(true),
      "Boolean: make EPS image of input and output?"
      )(
      "density_to_eps,d",
      value<bool>(&density_to_eps)
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

  CartogramInfo cart_info(visual_file_name,
                          world,
                          density_to_eps);

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

  // Determining name of input map
  std::string map_name = geo_file_name;
  if (map_name.find_last_of("/\\") != std::string::npos) {
    map_name = map_name.substr(map_name.find_last_of("/\\") + 1);
  }
  if (map_name.find('.') != std::string::npos) {
    map_name = map_name.substr(0, map_name.find('.'));
  }

  for (auto &inset_state : *cart_info.ref_to_inset_states()) {

    // Determining the name of the inset
    std::string inset_name = map_name;

    // Printing Inset Position if multiple insets present
    if (cart_info.n_insets() > 1) {
      inset_name = inset_name + "_" + inset_state.pos();
      std::cout << std::endl << std::endl
                << "Working on Inset with position: "
                << inset_state.pos()
                << std::endl;
    }
    inset_state.set_inset_name(inset_name);

    // Error checking Geometry
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

    // Rescale map to fit into a rectangular box [0, lx] * [0, ly].
    rescale_map(long_grid_side_length,
                &inset_state,
                cart_info.is_world_map());

    // Setting initial area errors
    inset_state.set_area_errs();

    // Filling density to fill horizontal adjacency map
    fill_with_density(&inset_state,
                      cart_info.trigger_write_density_to_eps());

    // Automatically coloring if no colors provided
    if (inset_state.colors_empty()) {
      auto_color(&inset_state);
    }

    // Writing EPS, if requested by command line option
    if (polygons_to_eps) {
      std::cout << "Writing " << inset_name << "_input.eps" << std::endl;
      write_map_to_eps((inset_name + "_input.eps"), &inset_state);
    }

    // Start map integration
    while (inset_state.n_finished_integrations() < max_integrations &&
           inset_state.max_area_err() > max_permitted_area_error) {

      std::cout << "Integration number "
                << inset_state.n_finished_integrations()
                << std::endl;


      if (inset_state.n_finished_integrations()  >  1) {
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
      inset_state.inc_integration();

      // Updating area errors
      inset_state.set_area_errs();
    }

    // Printing final cartogram
    json cart_json = cgal_to_json(&inset_state);
    write_to_json(cart_json,
                  geo_file_name,
                  (inset_name + "_cartogram_scaled.geojson"));

    // Printing EPS of output cartogram
    if (polygons_to_eps) {
      std::cout << "Writing " << inset_name << "_output.eps" << std::endl;
      write_map_to_eps((inset_name + "_output.eps"), &inset_state);
    }

    // Removing transformations
    unscale_map(&inset_state);

    // Printing unscaled cartogram
    cart_json = cgal_to_json(&inset_state);
    write_to_json(cart_json,
                  geo_file_name,
                  (inset_name + "_cartogram_unscaled.geojson"));

  }

  return EXIT_SUCCESS;
}
