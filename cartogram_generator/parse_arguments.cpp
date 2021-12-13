#include <iostream>
#include "constants.h"
#include "argparse.hpp"

argparse::ArgumentParser parsed_arguments(const int argc,
                                         const char *argv[],
                                         std::string &geo_file_name,
                                         std::string &visual_file_name,
                                         unsigned int &long_grid_side_length,
                                         bool &world,
                                         bool &triangulation,
                                         bool &make_csv,
                                         bool &make_polygon_eps,
                                         bool &output_equal_area,
                                         bool &output_to_stdout,
                                         bool &plot_density)
{

  // Creating parser for arguments using argparse.
  // From: https://github.com/p-ranav/argparse
  argparse::ArgumentParser arguments("./cartogram", "1.0");

  // Positional argument accepting geometry file (GeoJSON, JSON) as input.
  arguments.add_argument("geometry_file")
    .help("File path: GeoJSON file");

  // Optional argument accepting visual variables file (CSV) as input.
  arguments.add_argument("-v", "--visual_variable_file")
    .help("File path: CSV file with ID, area, and (optionally) colour");

  // Optional argument accepting long grid side length (unsigned int) as input.
  // Default value declared in "constants.h"
  arguments.add_argument("-l", "--long_grid_side_length")
    .default_value(default_long_grid_side_length)
    .scan<'u', unsigned int>()
    .help(
      "Integer: Number of grid cells along longer Cartesian coordinate axis"
    );

  // Optional boolean arguments.
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

  arguments.add_argument("-q", "--output_equal_area")
    .help("Boolean: Output equal area GeoJSON")
    .default_value(false)
    .implicit_value(true);

  arguments.add_argument("-t", "--triangulation")
    .help("Boolean: Project the cartogram using the triangulation?")
    .default_value(false)
    .implicit_value(true);

  arguments.add_argument("-m", "--make_csv")
    .help("Boolean: create CSV file from given GeoJSON?")
    .default_value(false)
    .implicit_value(true);

  arguments.add_argument("-s", "--output_to_stdout")
    .help("Boolean: Output GeoJSON to stdout")
    .default_value(false)
    .implicit_value(true);

  // Arguments regarding column names in provided visual variables file (CSV).
  arguments.add_argument("-i", "--id")
    .help(
      "String: Column name for IDs of geographic divisions [default: 1st CSV column]"
    );

  arguments.add_argument("-a", "--area")
    .help("String: Column name for target areas (default: 2nd CSV column)");

  arguments.add_argument("-c", "--color")
    .default_value(std::string("Color"))
    .help("String: Column name for colors");

  arguments.add_argument("-n", "--inset")
    .default_value(std::string("Inset"))
    .help("String: Column name for insets");

  // Parsing command line arguments.
  try {
    arguments.parse_args(argc, argv);
  }
  catch (const std::runtime_error& err) {
    std::cerr << "ERROR: " << err.what() << std::endl;
    std::cerr << arguments;
    std::exit(1);
  }

  // Setting long grid-side length.
  long_grid_side_length = arguments.get<unsigned int>("-l");

  // Setting boolean values.
  world = arguments.get<bool>("-w");
  triangulation = arguments.get<bool>("-t");
  make_csv = arguments.get<bool>("-m");
  make_polygon_eps = arguments.get<bool>("-e");
  output_equal_area = arguments.get<bool>("-q");
  output_to_stdout = arguments.get<bool>("-s");
  plot_density =  arguments.get<bool>("-d");

  // Printing geometry and visual file used.
  geo_file_name = arguments.get<std::string>("geometry_file");
  std::cerr << "Using geometry from file " << geo_file_name << std::endl;

  // Checking if visual variables file passed, or -m flag passed.
  if (arguments.present<std::string>("-v")) {
    visual_file_name = arguments.get<std::string>("-v");
    std::cerr << "Using visual variables from file "
              << visual_file_name
              << std::endl;
  } else if (!make_csv) {

    // CSV file not given, and user did not want to create one.
    std::cerr << "ERROR: No CSV file given!" << std::endl;
    std::cerr << "To create a CSV, please use the -m flag." << std::endl;
    std::cerr << arguments << std::endl;
    _Exit(15);
  }
  return arguments;
}