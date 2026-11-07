#include "parse_arguments.hpp"
#include "constants.hpp"
#include <cmath>
#include <limits>
#include <utility>

Arguments parse_arguments(const int argc, const char *argv[])
{
  Arguments args;

  // Create parser for arguments using argparse.
  // From https://github.com/p-ranav/argparse
  argparse::ArgumentParser arguments("./cartogram", version);

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
      "Integer: Number of starting grid cells along longer Cartesian "
      "coordinate axis");

  arguments.add_argument("-T", "--timeout")
    .default_value(std::numeric_limits<unsigned int>::max())
    .scan<'u', unsigned int>()
    .help(
      "Integer: Maximum time (in seconds) allowed for cartogram generation");

  arguments.add_argument("-N", "--max_allowed_autoscale_grid_length")
    .default_value(max_allowed_autoscale_grid_length)
    .scan<'u', unsigned int>()
    .help(
      "Integer: Maximum allowed number of grid cells along longer Cartesian "
      "coordinate axis");

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
    .help(
      "Boolean: Transform input GeoJSON into cartesian coordinates and exit")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("-S", "--disable_simplify_and_densify")
    .help(
      "Boolean: Disable iterative simplification and densification of "
      "polygons")
    .default_value(false)
    .implicit_value(true);
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

  arguments.add_argument("--do_not_fail_on_intersections")
    .help(
      "Boolean: Whether to still produce cartogram if polygons do not remain "
      "simple")
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
  arguments.add_argument("--id").help(
    pre + "IDs of geographic divisions [default: 1st CSV column header]");
  arguments.add_argument("--area").help(
    pre + "target areas [default: 2nd CSV column]");
  arguments.add_argument("--color", "--colour")
    .default_value(std::string("Color"))
    .help(pre + "colors");
  arguments.add_argument("--label")
    .default_value(std::string("Label"))
    .help(pre + "labels");
  arguments.add_argument("--inset")
    .default_value(std::string("Inset"))
    .help(pre + "insets");
  arguments.add_argument("--min_integrations")
    .help(
      "Integer: minimum number of integrations regardless of area error "
      "reached")
    .default_value(static_cast<unsigned int>(0))
    .scan<'u', unsigned int>();
  arguments.add_argument("--verbose")
    .help("Boolean: Print verbose time tracking output")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("--max_permitted_area_error")
    .help("Double: Maximum permitted area error for cartogram transformation")
    .default_value(default_max_permitted_area_error)
    .scan<'g', double>();
  arguments.add_argument("--quadtree_leaf_count_factor")
    .help("Unsigned int: Quadtree leaf count factor (should be a power of 2)")
    .default_value(default_quadtree_leaf_count_factor)
    .scan<'u', unsigned int>();
  arguments.add_argument("--conformality")
    .help(
      "Boolean: Track the cumulative map for the conformality postprocessor "
      "(experimental)")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("--conf_quad")
    .help(
      "Unsigned int: Conformality quadrature nodes per lattice cell along "
      "each axis")
    .default_value(static_cast<unsigned int>(2))
    .scan<'u', unsigned int>();
  arguments.add_argument("--conf_eval_quad")
    .help(
      "Unsigned int: Nodes per lattice cell along each axis of the "
      "independent conformality evaluation cloud")
    .default_value(static_cast<unsigned int>(4))
    .scan<'u', unsigned int>();
  arguments.add_argument("--conf_optimize")
    .help(
      "Boolean: Run the conformality optimizer on the tracked map "
      "(experimental; implies --conformality)")
    .default_value(false)
    .implicit_value(true);
  arguments.add_argument("--conf_modes")
    .help(
      "Unsigned int: Fourier modes of the conformality correction along the "
      "longer lattice axis")
    .default_value(static_cast<unsigned int>(16))
    .scan<'u', unsigned int>();
  arguments.add_argument("--conf_max_iter")
    .help(
      "Unsigned int: Maximum outer iterations of the conformality optimizer")
    .default_value(static_cast<unsigned int>(200))
    .scan<'u', unsigned int>();
  arguments.add_argument("--conf_alpha")
    .help(
      "Double: Exponent alpha in [0, 1] of the weight g0^alpha in the "
      "conformality objective, 0 for area on the input map, 1 for area on "
      "the cartogram")
    .default_value(0.0)
    .scan<'g', double>();
  arguments.add_argument("--conf_ocean_weight")
    .help(
      "Double: Factor (>= 0) on the conformality objective's weight of "
      "points outside every region")
    .default_value(1.0)
    .scan<'g', double>();
  arguments.add_argument("--conf_stall_tol")
    .help(
      "Double: The conformality optimizer stops when the last 10 accepted "
      "steps lowered the energy by at most this fraction (0: never)")
    .default_value(1e-5)
    .scan<'g', double>();
  arguments.add_argument("--conf_perturb")
    .help(
      "Double: Strain bound of a random smooth correction applied before the "
      "conformality optimization (0: none)")
    .default_value(0.0)
    .scan<'g', double>();
  arguments.add_argument("--conf_seed")
    .help("Unsigned int: Seed of the random correction (see --conf_perturb)")
    .default_value(static_cast<unsigned int>(1))
    .scan<'u', unsigned int>();

  // Parse command-line arguments
  try {
    arguments.parse_args(argc, argv);
  } catch (const std::runtime_error &err) {
    std::cerr << "ERROR: " << err.what() << ". ";
    std::cerr << arguments;
    std::exit(1);
  }

  // Ensure quadtree_leaf_count_factor is a power of 2
  unsigned int qlcf =
    arguments.get<unsigned int>("--quadtree_leaf_count_factor");
  // Round up to next power of 2 if not already a power of 2
  if ((qlcf <= 0) || ((qlcf & (~qlcf + 1)) != qlcf)) {
    std::cerr
      << "ERROR: --quadtree_leaf_count_factor must be an integer power of 2."
      << std::endl;
    std::exit(15);
  }
  args.quadtree_leaf_count_factor = qlcf;

  // Set long grid-side length
  args.n_grid_rows_or_cols = arguments.get<unsigned int>("-n");

  // Set maximum allowed grid length
  args.max_allowed_autoscale_grid_length =
    arguments.get<unsigned int>("--max_allowed_autoscale_grid_length");

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
  args.disable_simplification_densification =
    arguments.get<bool>("--disable_simplify_and_densify");
  args.remove_tiny_polygons = arguments.get<bool>("--remove_tiny_polygons");
  args.min_polygon_area = arguments.get<double>("--minimum_polygon_area");
  args.max_permitted_area_error =
    arguments.get<double>("--max_permitted_area_error");

  args.skip_projection = arguments.get<bool>("--skip_projection");
  args.make_csv = arguments.get<bool>("--make_csv");
  args.output_equal_area_map = arguments.get<bool>("--output_equal_area_map");
  args.redirect_exports_to_stdout =
    arguments.get<bool>("--redirect_exports_to_stdout");
  args.export_preprocessed = arguments.get<bool>("--export_preprocessed");
  args.export_time_report = arguments.get<bool>("--export_time_report");
  args.plot_density = arguments.get<bool>("--plot_density");
  args.plot_grid = arguments.get<bool>("--add_grid");
  args.plot_intersections = arguments.get<bool>("--plot_intersections");
  args.plot_polygons = arguments.get<bool>("--plot_polygons");
  args.plot_quadtree = arguments.get<bool>("--plot_quadtree");
  args.output_shifted_insets = arguments.get<bool>("--output_shifted_insets");
  args.do_not_fail_on_intersections =
    arguments.get<bool>("--do_not_fail_on_intersections");

  args.timeout_in_seconds = arguments.get<unsigned int>("--timeout");

  args.conf_optimize = arguments.get<bool>("--conf_optimize");
  args.conformality =
    arguments.get<bool>("--conformality") || args.conf_optimize;
  args.conf_quad = arguments.get<unsigned int>("--conf_quad");
  args.conf_eval_quad = arguments.get<unsigned int>("--conf_eval_quad");
  args.conf_modes = arguments.get<unsigned int>("--conf_modes");
  args.conf_max_iter = arguments.get<unsigned int>("--conf_max_iter");
  if (
    args.conformality && (args.conf_quad == 0 || args.conf_eval_quad == 0 ||
                          args.conf_modes == 0)) {
    std::cerr << "ERROR: --conf_quad, --conf_eval_quad and --conf_modes must "
                 "be positive."
              << std::endl;
    std::exit(17);
  }
  args.conf_alpha = arguments.get<double>("--conf_alpha");
  args.conf_ocean_weight = arguments.get<double>("--conf_ocean_weight");
  args.conf_stall_tol = arguments.get<double>("--conf_stall_tol");
  args.conf_perturb = arguments.get<double>("--conf_perturb");
  args.conf_seed = arguments.get<unsigned int>("--conf_seed");
  if (!(args.conf_alpha >= 0.0 && args.conf_alpha <= 1.0)) {
    std::cerr << "ERROR: --conf_alpha must be in [0, 1]." << std::endl;
    std::exit(17);
  }
  for (const auto &[flag, value] :
       {std::pair{"--conf_ocean_weight", args.conf_ocean_weight},
        std::pair{"--conf_stall_tol", args.conf_stall_tol},
        std::pair{"--conf_perturb", args.conf_perturb}}) {
    if (!(value >= 0.0) || !std::isfinite(value)) {
      std::cerr << "ERROR: " << flag << " must be finite and non-negative."
                << std::endl;
      std::exit(17);
    }
  }

  // arguments.present returns an optional
  args.id_col = arguments.present<std::string>("--id");
  args.area_col = arguments.present<std::string>("--area");
  args.inset_col = arguments.get<std::string>("--inset");
  args.color_col = arguments.get<std::string>("--color");
  args.label_col = arguments.get<std::string>("--label");

  // Check if user wants verbose output
  args.verbose = arguments.get<bool>("--verbose");

  // Check whether n_points is specified but --simplify_and_densify not passed
  if (args.disable_simplification_densification) {
    std::cerr << "WARNING: Simplification and densification disabled! "
              << "Polygons will not simplified (or densified). "
              << "This may result and in polygon intersections. "
              << "Thus, we are turning off topology checks. "
              << "To enable simplification, pass the -S flag." << std::endl;
    if (arguments.is_used("--n_points")) {
      std::cerr << "--n_points ignored." << std::endl;
    }
  }

  // Print names of geometry file
  if (arguments.is_used("geometry_file")) {
    args.geo_file_name = arguments.get<std::string>("geometry_file");
    std::cerr << "Using geometry from file " << args.geo_file_name
              << std::endl;
  } else {

    // GeoJSON file not provided
    std::cerr << arguments << std::endl;
    std::cerr << "ERROR: No Geometry file provided!" << std::endl;
    std::cerr << "Please provide a geometry file in standard GeoJSON format."
              << std::endl;
    std::exit(16);
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
    std::exit(15);
  } else {
    args.visual_file_name = "";
  }

  return args;
}