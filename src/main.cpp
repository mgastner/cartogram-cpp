#include "cartogram_info.hpp"
#include "parse_arguments.hpp"
#include "progress_tracker.hpp"

int main(const int argc, const char *argv[])
{
  // Struct to store command line arguments
  Arguments args;

  // Parse command-line arguments
  argparse::ArgumentParser arguments = parsed_arguments(argc, argv, args);

  // Initialize cart_info. It contains all the information about the cartogram
  // that needs to be handled by functions called from main().
  CartogramInfo cart_info(args);

  // Read geometry. If the GeoJSON does not explicitly contain a "crs" field,
  // we assume that the coordinates are in longitude and latitude.
  cart_info.read_geojson();

  if (!args.visual_file_name.empty()) {

    // Read visual variables (e.g., area and color) from CSV
    cart_info.read_csv(arguments);
  }

  // Project to equal area, if necessary
  // Write Output and EXIT if --output_equal_area_map is set to true
  cart_info.project_to_equal_area();

  // Store total number of GeoDivs to monitor progress
  double total_geo_divs = cart_info.n_geo_divs();

  // Write input map, with insets nicely placed
  if (args.plot_polygons) {
    cart_info.plot_input();
  }

  // Move insets to desired positions and EXIT
  if (args.output_shifted_insets) {
    cart_info.write_shifted_insets();
  }

  // Track progress of the cartogram generation
  ProgressTracker progress_tracker(total_geo_divs);

  // Preprocess Insets for Integration:
  // -- Set inset name: map_name + "_" + inset_pos
  // -- Rescale
  // -- Replace missing and zero target areas
  // -- Simplify
  // -- Remove tiny polygons
  // -- Write input map if requested (and color polygons, if necessary)
  cart_info.preprocess();

  // Iterate over insets and integrate each one
  for (InsetState &inset_state : cart_info.ref_to_inset_states()) {
    inset_state.integrate(progress_tracker);
  }

  // Rescale insets in correct proportion to each other
  cart_info.rescale_insets();

  // Shift insets so that they do not overlap
  cart_info.reposition_insets(args.redirect_exports_to_stdout);

  // Output to GeoJSON
  cart_info.write_geojson("cartogram");

  if (args.plot_polygons) {
    cart_info.write_svg("cartogram");
  }

  // Stop total time timer, and print time summary report
  // Export report with area errors to CSV if requested
  cart_info.print_time_report();
  return EXIT_SUCCESS;
}