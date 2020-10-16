// TO DO: positional matching of argument flags

#include "constants.h"
#include "map_state.h"
#include "read_csv.h"
#include "read_geojson.h"
#include "rescale_map.h"
#include <boost/program_options.hpp>
#include <iostream>

// Functions that are called if the corresponding command-line options are
// present
void on_geometry(const std::string geometry_file_name)
{
  std::cerr << "Using geometry from file " << geometry_file_name << std::endl;
  return;
}

void on_visual_variables(const std::string geometry_file_name)
{
  std::cerr << "Using visual variables from file "
            << geometry_file_name
            << std::endl;
  return;
}

int main(const int argc, const char *argv[])
{
  using namespace boost::program_options;
  std::string geo_file_name;

  // Default number of grid cells along longer Cartesian coordinate axis.
  int long_lattice_side_length = default_long_lattice_side_length;

  // World maps need special projections. By default, we assume that the
  // input map is not a world map.
  bool world = false;

  // Parse command-line options. See
  // https://theboostcpplibraries.com/boost.program_options
  try {
    options_description desc{"Options"};
    desc.add_options()(
      "help,h", "Help screen"
      )(
      "geometry,g",
      value<std::string>(&geo_file_name)->required()->notifier(on_geometry),
      "GeoJSON file"
      )(
      "visual_variables,v",
      value<std::string>()->notifier(on_visual_variables),
      "CSV file with area and (optionally) colour"
      )(
      "long_lattice_side_length,l",
      value<int>(&long_lattice_side_length),
      "Number of grid cells along longer Cartesian coordinate axis"
      )(
      "world,w",
      value<bool>(&world),
      "Boolean: is input a world map in longitude-latitude format?"
      );
    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    if (vm.count("help") || vm.empty()) {
      std::cerr << desc << '\n';
      return EXIT_SUCCESS;
    } else {
      notify(vm);  // Triggers notifier functions such as on_geometry()
    }
  } catch (const error &ex) {
    std::cerr << "ERROR: " << ex.what() << std::endl;
    return EXIT_FAILURE;
  }
  if ((long_lattice_side_length <= 0) ||
      ((long_lattice_side_length &
        (~long_lattice_side_length + 1)) != long_lattice_side_length)) {
    std::cerr << "ERROR: long_lattice_side_length must be an integer"
              << "power of 2." << std::endl;
    _Exit(15);
  }
  MapState map_state(world);

  // Read visual variables (e.g. area) from CSV
  read_csv("What's", "up", "doc?");

  // Read geometry
  try {
    read_geojson(geo_file_name, &map_state);
  } catch (const std::system_error& e) {
    std::cerr << "ERROR: " << e.what() << " (" << e.code() << ")" << std::endl;
    return EXIT_FAILURE;
  }

  // Rescale map to fit into a rectangular box [0, lx] * [0, ly].
  rescale_map(long_lattice_side_length, &map_state);

  for (auto gd : map_state.get_geo_divs()) {
    for (auto pwh : gd.get_polygons_with_holes()) {
      CGAL::set_pretty_mode(std::cout);
      std::cout << pwh << std::endl;
    }
  }
  return EXIT_SUCCESS;
}
