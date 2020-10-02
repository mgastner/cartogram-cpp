// TO DO: positional matching of argument flags

#include "read_csv.h"
#include "read_geojson.h"
#include <boost/program_options.hpp>
#include <iostream>

// Functions that are called if the corresponding command-line options are
// present
void on_geometry(std::string geometry_file_name)
{
  std::cerr << "Using geometry from file " << geometry_file_name << std::endl;
  return;
}
void on_visual_variables(std::string geometry_file_name)
{
  std::cerr << "Using visual variables from file "
            << geometry_file_name
            << std::endl;
  return;
}

int main(int argc, const char *argv[])
{
  using namespace boost::program_options;

  std::string geo_file_name;

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
      "CSV file with area and (optionally) colour");
    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    if (vm.count("help") || vm.empty()) {
      std::cerr << desc << '\n';
      return EXIT_SUCCESS;
    } else {
      notify(vm);  // Triggers notifier functions such as on_geometry()
    }
  } catch (const error &ex) {
    std::cerr << "ERROR: " << ex.what() << '\n';
    return EXIT_FAILURE;
  }

  // Read visual variables (e.g. area) from CSV
  read_csv("What's", "up", "doc?");

  // Read geometry
  try {
    read_geojson(geo_file_name);
  } catch (const std::system_error& e) {
    std::cerr << "ERROR: " << e.what() << " (" << e.code() << ")" << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
