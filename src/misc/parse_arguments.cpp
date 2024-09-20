#include "parse_arguments.hpp"
#include "constants.hpp"

argparse::ArgumentParser parsed_arguments(
  const int argc,
  const char *argv[],
  std::string &geo_file_name,
  std::string &visual_file_name,
  unsigned int &long_grid_side_length,
  unsigned int &target_points_per_inset,
  bool &world,
  bool &triangulation,
  bool &qtdt_method,
  bool &simplify,
  bool &make_csv,
  bool &output_equal_area,
  bool &output_to_stdout,
  bool &plot_density,
  bool &plot_grid,
  bool &plot_intersections,
  bool &plot_polygons,
  bool &remove_tiny_polygons,
  double &minimum_polygon_area,
  bool &plot_quadtree,
  bool &rays)
{
  // Create parser for arguments using argparse.
  // From https://github.com/p-ranav/argparse
  argparse::ArgumentParser arguments("./cartogram", "2.0");

  // Positional argument accepting geometry file (GeoJSON, JSON) as input
  arguments.add_argument("geometry_file")
    .default_value("none")
    .help("File path: GeoJSON file");

  // Positional argument accepting visual variables file (CSV) as input
  arguments.add_argument("visual_variable_file")
    .default_value("none")
    .help("File path: CSV file with ID, area, and (optionally) colour");

  // Optional argument accepting long grid side length (unsigned int) as
  // input. Default value declared in "constants.hpp"
  arguments.add_argument("-n", "--n_grid_rows_or_cols")
    .default_value(default_long_grid_length)
    .scan<'u', unsigned int>()
    .help(
      "Integer: Number of grid cells along longer Cartesian coordinate axis");

  // Optional boolean arguments
  arguments.add_argument("-W", "--world")
    .help("Boolean: is input a world map in longitude-latitude format?")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-p", "--plot_polygons")
    .help("Boolean: Plot images of input and output cartogram")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-q", "--plot_quadtree")
    .help("Boolean: Plot images of Quadtree-Delaunay Triangulation")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-d", "--plot_density")
    .help("Boolean: Plot images of flatten and blur density")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-g", "--plot_grid")
    .help("Boolean: Plot images of with transformed grid grid")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-i", "--plot_intersections")
    .help("Boolean: Plot images of intersections (if any)")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-E", "--output_equal_area")
    .help(
      "Boolean: Output equal area GeoJSON from input GeoJSON - no cartogram")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-T", "--triangulation")
    .help("Boolean: Project the cartogram using the triangulation method")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-Q", "--qtdt_method")
    .help("Boolean: Use Quadtree-Delaunay Triangulation Method")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-S", "--simplify")
    .help("Boolean: Simplify polygons")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-P", "--n_points")
    .help(
      "Integer: If simplification enabled, target number of points per inset")
    .default_value(default_target_points_per_inset)
    .scan<'u', unsigned int>();
  arguments.add_argument("-M", "--make_csv")
    .help("Boolean: create CSV file from given GeoJSON?")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-O", "--output_to_stdout")
    .help("Boolean: Output GeoJSON to stdout")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-R", "--remove_tiny_polygons")
    .help("Boolean: Remove tiny polygons")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-m", "--minimum_polygon_size")
    .help(
      std::string("Double: If remove-tiny-polygons enabled, ") +
      "minimum size of polygons as proportion of total area")
    .default_value(default_minimum_polygon_area)
    .scan<'g', double>();
  arguments.add_argument("-r", "--use_ray_shooting_method")
    .help("Boolean: Use old ray shooting method to fill density")
    .default_value(false)
    .implicit_value(true);

  // Arguments of column names in provided visual variables file (CSV)
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

  // Parse command-line arguments
  try {
    arguments.parse_args(argc, argv);
  } catch (const std::runtime_error &err) {
    std::cerr << "ERROR: " << err.what() << std::endl;
    std::cerr << arguments;
    std::exit(1);
  }

  // Set long grid-side length
  long_grid_side_length = arguments.get<unsigned int>("-n");

  // If world flag is set, and long-gride side length is not explicitly set,
  // then 512 makes the output look better
  if (arguments.get<bool>("-W") && !arguments.is_used("-n")) {
    long_grid_side_length = 512;
  }

  // Set target_points_per_inset
  target_points_per_inset = arguments.get<unsigned int>("-P");

  // Set boolean values
  world = arguments.get<bool>("-W");
  triangulation = arguments.get<bool>("-T");
  qtdt_method = arguments.get<bool>("-Q");
  simplify = arguments.get<bool>("-S");
  remove_tiny_polygons = arguments.get<bool>("-R");
  minimum_polygon_area = arguments.get<double>("-m");
  rays = arguments.get<bool>("-r");
  if (!triangulation && simplify) {

    // If tracer points are on the FTReal2d, then simplification requires
    // triangulation. Otherwise, the cartogram may contain intersecting lines
    // before simplification. The intersection coordinates would not be
    // explicit in the non-simplified polygons, but they would be added to the
    // simplified polygon, making it more difficult to uniquely match
    // non-simplified and simplified polygons.
    triangulation = true;
  }
  make_csv = arguments.get<bool>("-M");
  output_equal_area = arguments.get<bool>("-E");
  output_to_stdout = arguments.get<bool>("-O");
  plot_density = arguments.get<bool>("-d");
  plot_grid = arguments.get<bool>("-g");
  plot_intersections = arguments.get<bool>("-i");
  plot_polygons = arguments.get<bool>("-p");
  plot_quadtree = arguments.get<bool>("-q");

  if (
    arguments.is_used("-O") && !arguments.is_used("-S") &&
    !arguments.is_used("-Q")) {
    std::cerr << "ERROR: --simplify flag not passed!\n";
    std::cerr << "--output_to_stdout flag is only supported with "
                 "simplification or quadtree.\n";
    std::cerr << "To enable simplification, pass the -S flag.\n";
    std::cerr << "To enable quadtree, pass the -Q flag.\n";
    std::cerr << arguments << std::endl;
    _Exit(18);
  }

  // Check whether n_points is specified but --simplify not passed
  if (arguments.is_used("-P") && !arguments.is_used("-S")) {
    std::cerr << "WARNING: --simplify flag not passed!" << std::endl;
    std::cerr << "Polygons will not be simplified." << std::endl;
    std::cerr << "To enable simplification, pass the -S flag." << std::endl;
    std::cerr << arguments << std::endl;
  }

  // Check whether T flag is set, but not Q
  if (arguments.is_used("-T") && !arguments.is_used("-Q")) {
    std::cerr << "ERROR: --qtdt_method flag not passed!" << std::endl;
    std::cerr << "QTDT method is necessary for Quadtree images." << std::endl;
    std::cerr << "To use qtdt method, pass the -Q flag." << std::endl;
    std::cerr << arguments << std::endl;
    _Exit(17);
  }

  // Print names of geometry file
  if (arguments.is_used("geometry_file")) {
    geo_file_name = arguments.get<std::string>("geometry_file");
    std::cerr << "Using geometry from file " << geo_file_name << std::endl;
  } else {

    // GeoJSON file not provided
    std::cerr << arguments << std::endl;
    std::cerr << "ERROR: No Geometry file provided!" << std::endl;
    std::cerr << "Please provide a geometry file in standard GeoJSON format."
              << std::endl;
    _Exit(16);
  }

  // Check if a visual-variables file or -m flag is passed
  if (arguments.is_used("visual_variable_file")) {
    visual_file_name = arguments.get<std::string>("visual_variable_file");
    std::cerr << "Using visual variables from file " << visual_file_name
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
