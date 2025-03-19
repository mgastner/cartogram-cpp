#include "parse_arguments.hpp"
#include "constants.hpp"

void parsed_arguments(
  const int argc,
  const char *argv[],
  Arguments &args)
{
  // Create parser for arguments using argparse.
  // From https://github.com/p-ranav/argparse
  argparse::ArgumentParser arguments("./cartogram", RELEASE_TAG);

  // Positional argument accepting geometry file (GeoJSON, JSON) as input
  arguments.add_argument("geometry_file")
    .required()
    .help("File path: GeoJSON file");

  // Positional argument accepting visual variables file (CSV) as input
  arguments.add_argument("visual_variable_file")
    .default_value(std::string("none"))
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
  arguments.add_argument("-g", "--add_grid")
    .help("Boolean: Add area legend and grid to relevant plots")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-i", "--plot_intersections")
    .help("Boolean: Plot images of intersections (if any)")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-E", "--output_equal_area_map")
    .help("Boolean: Transform input GeoJSON into cartesian coordinates and exit")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-T", "--triangulation")
    .help("Boolean: Enable cartogram projection via triangulation")
    .default_value(true)
    .implicit_value(false);
  arguments.add_argument("-Q", "--qtdt_method")
    .help("Boolean: Enable Quadtree-Delaunay Triangulation Method")
    .default_value(true)
    .implicit_value(false);
  arguments.add_argument("-S", "--simplify_and_densify")
    .help("Boolean: Enable iterative simplification and densification of polygons")
    .default_value(true)
    .implicit_value(false);
  arguments.add_argument("--skip_projection")
    .help("Boolean: Skip projection to equal area")
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
  arguments.add_argument("-O", "--redirect_exports_to_stdout")
    .help("Boolean: Redirect all exports to stdout as valid JSON")
    // Currently, this help message is not accurate.
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("--output_shifted_insets")
    .help(
      "Boolean: Output repositioned insets in cartesian coordinates GeoJSON")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-R", "--remove_tiny_polygons")
    .help("Boolean: Remove tiny polygons")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-m", "--minimum_polygon_area")
    .help("Double: Minimum size of tiny polygons as proportion of total area")
    .default_value(default_minimum_polygon_area)
    .scan<'g', double>();
  arguments.add_argument("-r", "--use_ray_shooting_method")
    .help("Boolean: Use old ray shooting method to fill density")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("--export_preprocessed")
    .help("Boolean: write input GeoJSON and CSV after preprocessing")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("--export_time_report")
    .help("Boolean: write extended time report to CSV file")
    .default_value(false)
    .implicit_value(true);

  // Arguments of column names in provided visual variables file (CSV)
  std::string pre = "String: Column name for ";
  arguments.add_argument("-D", "--id")
    .help(pre + "IDs of geographic divisions [default: 1st CSV column header]");
  arguments.add_argument("-A", "--area")
    .help(pre + "target areas [default: 2nd CSV column]");
  arguments.add_argument("-C", "--color", "--colour")
    .default_value(std::string("Color"))
    .help(pre + "colors");
  arguments.add_argument("-L", "--label")
    .default_value(std::string("Label"))
    .help(pre + "labels");
  arguments.add_argument("-I", "--inset")
    .default_value(std::string("Inset"))
    .help(pre + "insets");
  arguments.add_argument("--min_integrations")
    .help("Integer: minimum number of integrations regardless of area error reached")
    .default_value(static_cast<unsigned int>(0))
    .scan<'u', unsigned int>();

  // Parse command-line arguments
  try {
    arguments.parse_args(argc, argv);
  } catch (const std::runtime_error &err) {
    std::cerr << "ERROR: " << err.what() << ". ";
    std::cerr << arguments;
    std::exit(1);
  }

  // Set long grid-side length
  args.n_grid_rows_or_cols = arguments.get<unsigned int>("-n");

  // If world flag is set, and long-gride side length is not explicitly set,
  // then 512 makes the output look better
  if (arguments.get<bool>("--world") && !arguments.is_used("-n")) {
    args.n_grid_rows_or_cols = 512;
  }

  // Set target_points_per_inset
  args.target_points_per_inset = arguments.get<unsigned int>("-P");

  args.min_integrations = arguments.get<unsigned int>("--min_integrations");

  // Set boolean values
  args.world = arguments.get<bool>("--world");
  args.triangulation = arguments.get<bool>("--triangulation");
  args.qtdt_method = arguments.get<bool>("--qtdt_method");
  args.simplify = arguments.get<bool>("--simplify_and_densify");
  args.remove_tiny_polygons = arguments.get<bool>("--remove_tiny_polygons");
  args.min_polygon_area = arguments.get<double>("--minimum_polygon_area");
  args.rays = arguments.get<bool>("--use_ray_shooting_method");
  if (!args.triangulation && args.simplify) {

    // If tracer points are on the FTReal2d, then simplification requires
    // triangulation. Otherwise, the cartogram may contain intersecting lines
    // before simplification. The intersection coordinates would not be
    // explicit in the non-simplified polygons, but they would be added to the
    // simplified polygon, making it more difficult to uniquely match
    // non-simplified and simplified polygons.
    args.triangulation = true;
  }
  args.skip_projection = arguments.get<bool>("--skip_projection");
  args.make_csv = arguments.get<bool>("--make_csv");
  args.output_equal_area_map = arguments.get<bool>("--output_equal_area_map");
  args.redirect_exports_to_stdout = arguments.get<bool>("--redirect_exports_to_stdout");
  args.export_preprocessed = arguments.get<bool>("--export_preprocessed");
  args.export_time_report = arguments.get<bool>("--export_time_report");
  args.plot_density = arguments.get<bool>("--plot_density");
  args.plot_grid = arguments.get<bool>("--add_grid");
  args.plot_intersections = arguments.get<bool>("--plot_intersections");
  args.plot_polygons = arguments.get<bool>("--plot_polygons");
  args.plot_quadtree = arguments.get<bool>("--plot_quadtree");
  args.output_shifted_insets =
    arguments.get<bool>("--output_shifted_insets");

  args.id_col = arguments.present<std::string>("--id");
  args.area_col = arguments.present<std::string>("--area");
  args.inset_col = arguments.get<std::string>("--inset");
  args.color_col = arguments.get<std::string>("--color");
  args.label_col = arguments.get<std::string>("--label");
  // Check if user wants to redirect output to stdout
  if (arguments.is_used("-O") && !args.simplify && !args.qtdt_method) {
    std::cerr << "ERROR: simplification disabled!\n";
    std::cerr << "--output_to_stdout flag is only supported with "
                 "simplification or quadtree.\n";
    std::cerr << "To enable simplification, do not pass the -S flag.\n";
    std::cerr << "To enable quadtree, do not pass the -Q flag.\n";
    std::cerr << arguments << std::endl;
    _Exit(18);
  }

  // Check whether n_points is specified but --simplify_and_densify not passed
  if (!args.simplify) {
    std::cerr << "WARNING: Simplification and densification disabled! "
              << "Polygons will not simplified (or densified). "
              << "This may result and in polygon intersections. "
              << "Thus, we are turning off topology checks. "
              << "To enable simplification, pass the -S flag." << std::endl;
    if (arguments.is_used("--n_points")) {
      std::cerr << "--n_points ignored." << std::endl;
    }
    std::cerr << arguments << std::endl;
  }

  // Check whether T flag is set, but not Q
  if (args.triangulation && !args.qtdt_method) {
    std::cerr
      << "ERROR: Can't disable qtdt_method without disabling triangulation."
      << std::endl;
    std::cerr << "QTDT method is necessary for Quadtree images." << std::endl;
    std::cerr << "To disable Triangulation, pass the -T flag." << std::endl;
    std::cerr << arguments << std::endl;
    _Exit(17);
  }

  // Print names of geometry file
  if (arguments.is_used("geometry_file")) {
    args.geo_file_name = arguments.get<std::string>("geometry_file");
    std::cerr << "Using geometry from file " << args.geo_file_name << std::endl;
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
    args.visual_file_name = arguments.get<std::string>("visual_variable_file");
    std::cerr << "Using visual variables from file " << args.visual_file_name
              << std::endl;
  } else if (!args.make_csv and !args.output_equal_area_map) {

    // CSV file not given, and user does not want to create one
    std::cerr << arguments << std::endl;
    std::cerr << "ERROR: No CSV file provided!" << std::endl;
    std::cerr << "To create a CSV, please use the -m flag." << std::endl;
    _Exit(15);
  } else {
    args.visual_file_name = "";
  }
}