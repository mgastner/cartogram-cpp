#include "constants.h"
#include "argparse.hpp"
#include <iostream>

argparse::ArgumentParser parsed_arguments(const int argc,
                                          const char *argv[],
                                          std::string &geo_file_name,
                                          std::string &visual_file_name,
                                          unsigned int &long_grid_side_length,
                                          unsigned int &target_points_per_inset,
                                          bool &world,
                                          bool &triangulation,
                                          bool &simplify,
                                          bool &make_csv,
                                          bool &make_polygon_eps,
                                          bool &output_equal_area,
                                          bool &output_to_stdout,
                                          bool &plot_density,
                                          bool &plot_graticule,
                                          bool &plot_intersections)
{
  // Create parser for arguments using argparse.
  // From https://github.com/p-ranav/argparse
  argparse::ArgumentParser arguments("./cartogram", "1.0");

  // Positional argument accepting geometry file (GeoJSON, JSON) as input
  arguments.add_argument("geometry_file")
  .help("File path: GeoJSON file");

  // Positional argument accepting visual variables file (CSV) as input
  arguments.add_argument("visual_variable_file")
  .default_value("none")
  .help("File path: CSV file with ID, area, and (optionally) colour");

  // Optional argument accepting long grid side length (unsigned int) as
  // input. Default value declared in "constants.h"
  arguments.add_argument("-l", "--long_grid_side_length")
  .default_value(default_long_grid_side_length)
  .scan<'u', unsigned int>()
  .help(
    "Integer: Number of grid cells along longer Cartesian coordinate axis");

  // Optional boolean arguments
  arguments.add_argument("-w", "--world")
  .help("Boolean: is input a world map in longitude-latitude format?")
  .default_value(false)
  .implicit_value(true);

  arguments.add_argument("-e", "--polygons_to_eps")
  .help("Boolean: make EPS image of input and output?")
  .default_value(false)
  .implicit_value(true);

  arguments.add_argument("-d", "--density_to_eps")
  .help("Boolean: make EPS images *_density_*.eps?")
  .default_value(false)
  .implicit_value(true);

  arguments.add_argument("-g", "--graticule_to_eps")
  .help("Boolean: make EPS images with graticule?")
  .default_value(false)
  .implicit_value(true);

  arguments.add_argument("-I", "--intersections_to_eps")
  .help("Boolean: make EPS images *_intersections_*.eps?")
  .default_value(false)
  .implicit_value(true);

  arguments.add_argument("-q", "--output_equal_area")
  .help("Boolean: Output equal area GeoJSON")
  .default_value(false)
  .implicit_value(true);

  arguments.add_argument("-t", "--triangulation")
  .help("Boolean: Project the cartogram using the triangulation method?")
  .default_value(false)
  .implicit_value(true);

  arguments.add_argument("-s", "--simplify")
  .help("Boolean: Shall the polygons be simplified?")
  .default_value(false)
  .implicit_value(true);

  arguments.add_argument("-S", "--n_points")
  .help("Integer: If simplification enabled, target number of points per inset")
  .default_value(default_target_points_per_inset)
  .scan<'u', unsigned int>();

  arguments.add_argument("-m", "--make_csv")
  .help("Boolean: create CSV file from given GeoJSON?")
  .default_value(false)
  .implicit_value(true);

  arguments.add_argument("-o", "--output_to_stdout")
  .help("Boolean: Output GeoJSON to stdout")
  .default_value(false)
  .implicit_value(true);

  // Arguments regarding column names in provided visual variables file (CSV)
  std::string pre = "String: Column name for ";
  arguments.add_argument("-i", "--id")
  .help(pre + "IDs of geographic divisions [default: 1st CSV column]");

  arguments.add_argument("-a", "--area")
  .help(pre + "target areas [default: 2nd CSV column]");

  arguments.add_argument("-c", "--color")
  .default_value(std::string("Color"))
  .help(pre + "colors");

  arguments.add_argument("-n", "--inset")
  .default_value(std::string("Inset"))
  .help(pre + "insets");

  // Parse command line arguments
  try {
    arguments.parse_args(argc, argv);
  }
  catch (const std::runtime_error& err) {
    std::cerr << "ERROR: " << err.what() << std::endl;
    std::cerr << arguments;
    std::exit(1);
  }

  // Set long grid-side length
  long_grid_side_length = arguments.get<unsigned int>("-l");

  // Set target_points_per_inset
  target_points_per_inset = arguments.get<unsigned int>("-S");

  // Set boolean values
  world = arguments.get<bool>("-w");
  triangulation = arguments.get<bool>("-t");
  simplify = arguments.get<bool>("-s");
  if (!triangulation && simplify) {

    // Simplification requires triangulation. Otherwise, the cartogram may
    // contain intersecting lines before simplification. The intersection
    // coordinates would not be explicit in the non-simplified polygons,
    // but they would be added to the simplified polygon, making it more
    // difficult to uniquely match non-simplified and simplified polygons.
    // triangulation = true;
  }
  make_csv = arguments.get<bool>("-m");
  make_polygon_eps = arguments.get<bool>("-e");
  output_equal_area = arguments.get<bool>("-q");
  output_to_stdout = arguments.get<bool>("-o");
  plot_density =  arguments.get<bool>("-d");
  plot_graticule = arguments.get<bool>("-g");
  plot_intersections = arguments.get<bool>("-I");

  // Check if n_polygons specified, but --simplify not passed.
  if (arguments.is_used("-S") && !arguments.is_used("-s")) {
    std::cerr << "WARNING: --simplify flag not passed!" << std::endl;
    std::cerr << "Polygons will not be simplified." << std::endl;
    std::cerr << "To enable simplification, pass the -s flag." << std::endl;
    std::cerr << arguments << std::endl;
  }

  // Print geometry and visual-variables file used
  geo_file_name = arguments.get<std::string>("geometry_file");
  std::cerr << "Using geometry from file " << geo_file_name << std::endl;

  // Check if a visual-variables file or -m flag is passed
  if (arguments.is_used("visual_variable_file")) {
    visual_file_name = arguments.get<std::string>("visual_variable_file");
    std::cerr << "Using visual variables from file "
              << visual_file_name
              << std::endl;
  } else if (!make_csv) {

    // CSV file not given, and user does not want to create one
    std::cerr << arguments << std::endl;
    std::cerr << "ERROR: No CSV file given!" << std::endl;
    std::cerr << "To create a CSV, please use the -m flag." << std::endl;
    _Exit(15);
  }
  return arguments;
}
