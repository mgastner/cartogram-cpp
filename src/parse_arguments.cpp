#include "constants.h"
#include "argparse.hpp"
#include <iostream>

argparse::ArgumentParser parsed_arguments(
  const int argc,
  const char *argv[],
  std::string &geo_file_name,
  std::string &visual_file_name,
  unsigned int
  &max_n_graticule_rows_or_cols,
  unsigned int &target_points_per_inset,
  bool &world,
  bool &triangulation,
  bool &simplify,
  bool &make_csv,
  bool &produce_map_image,
  bool &image_format_ps,
  bool &output_equal_area,
  bool &output_to_stdout,
  bool &plot_density,
  bool &plot_graticule,
  bool &plot_graticule_heatmap,
  bool &plot_intersections)
{
  // Create parser for arguments using argparse.
  // From https://github.com/p-ranav/argparse
  argparse::ArgumentParser arguments("./cartogram", "1.0");

  // Positional argument accepting geometry file (GeoJSON, JSON) as input
  arguments.add_argument("geometry_file")
  .default_value("none")
  .help("File path: GeoJSON file");

  // Positional argument accepting visual variables file (CSV) as input
  arguments.add_argument("visual_variable_file")
  .default_value("none")
  .help("File path: CSV file with ID, area, and (optionally) colour");

  // Optional argument accepting long grid side length (unsigned int) as
  // input. Default value declared in "constants.h"
  arguments.add_argument("-N", "--n_graticule_rows_or_cols")
  .default_value(default_max_n_graticule_rows_or_cols)
  .scan<'u', unsigned int>()
  .help(
    "Integer: Number of grid cells along longer Cartesian coordinate axis");

  // Optional boolean arguments
  arguments.add_argument("-w", "--world")
  .help("Boolean: is input a world map in longitude-latitude format?")
  .default_value(false)
  .implicit_value(true);

  arguments.add_argument("-e", "--map_image")
  .help("Boolean: produce SVG/PS image of input and output?")
  .default_value(false)
  .implicit_value(true);
  
  arguments.add_argument("-p", "--image_format_ps")
  .help("Boolean: use .ps format for images?")
  .default_value(false)
  .implicit_value(true);

  arguments.add_argument("-d", "--density_image")
  .help("Boolean: produce density images *_density_*.svg/ps?")
  .default_value(false)
  .implicit_value(true);

  arguments.add_argument("-g", "--add_graticules_to_image")
  .help("Boolean: include graticules in images?")
  .default_value(false)
  .implicit_value(true);
  
  arguments.add_argument("-h", "--graticule_heatmap_image")
  .help("Boolean: produce graticule heatmap images *_graticule_heatmap.svg/ps?")
  .default_value(false)
  .implicit_value(true);

  arguments.add_argument("-i", "--intersections_image")
  .help("Boolean: produce intersections images *_intersections_*.svg/ps?")
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

  arguments.add_argument("-P", "--n_points")
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
  arguments.add_argument("-D", "--id")
  .help(pre + "IDs of geographic divisions [default: 1st CSV column]");

  arguments.add_argument("-A", "--area")
  .help(pre + "target areas [default: 2nd CSV column]");

  arguments.add_argument("-C", "--color")
  .default_value(std::string("Color"))
  .help(pre + "colors");

  arguments.add_argument("-L", "--label")
  .default_value(std::string("Label"))
  .help(pre + "labels");

  arguments.add_argument("-I", "--inset")
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
  max_n_graticule_rows_or_cols = arguments.get<unsigned int>("-N");

  // Set target_points_per_inset
  target_points_per_inset = arguments.get<unsigned int>("-P");

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
    triangulation = true;
  }
  make_csv = arguments.get<bool>("-m");
  produce_map_image = arguments.get<bool>("-e");
  image_format_ps = arguments.get<bool>("-p");
  output_equal_area = arguments.get<bool>("-q");
  output_to_stdout = arguments.get<bool>("-o");
  plot_density =  arguments.get<bool>("-d");
  plot_graticule = arguments.get<bool>("-g");
  plot_graticule_heatmap = arguments.get<bool>("-h");
  plot_intersections = arguments.get<bool>("-i");

  // Check whether n_points is specified but --simplify not passed
  if (arguments.is_used("-P") && !arguments.is_used("-s")) {
    std::cerr << "WARNING: --simplify flag not passed!" << std::endl;
    std::cerr << "Polygons will not be simplified." << std::endl;
    std::cerr << "To enable simplification, pass the -s flag." << std::endl;
    std::cerr << arguments << std::endl;
  }

  // Print names of geometry and visual-variables files used
  if (arguments.is_used("geometry_file")) {
    geo_file_name = arguments.get<std::string>("geometry_file");
    std::cerr << "Using geometry from file " << geo_file_name << std::endl;

  } else {

    // GeoJSON file not provided
    std::cerr << arguments << std::endl;
    std::cerr << "ERROR: No Geometry file provided!" << std::endl;
    std::cerr <<
      "Please provide a geometry file in the standard GeoJSON format."
              << std::endl;
    _Exit(16);
  }

  // Check if a visual-variables file or -m flag is passed
  if (arguments.is_used("visual_variable_file")) {
    visual_file_name = arguments.get<std::string>("visual_variable_file");
    std::cerr << "Using visual variables from file "
              << visual_file_name
              << std::endl;
  } else if (!make_csv) {

    // CSV file not given, and user does not want to create one
    std::cerr << arguments << std::endl;
    std::cerr << "ERROR: No CSV file provided!" << std::endl;
    std::cerr << "To create a CSV, please use the -m flag." << std::endl;
    _Exit(15);
  }
  return arguments;
}
