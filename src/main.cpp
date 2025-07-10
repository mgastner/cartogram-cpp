#include "cartogram_info.hpp"
#include "parse_arguments.hpp"
#include "progress_tracker.hpp"

int main(const int argc, const char *argv[])
{
  // Parse command-line arguments
  Arguments args = parse_arguments(argc, argv);

  // Initialize cart_info. It contains all the information about the cartogram
  // that needs to be handled by functions called from main().
  CartogramInfo cart_info(args);

  // Read geometry. If the GeoJSON does not explicitly contain a "crs" field,
  // we assume that the coordinates are in longitude and latitude.
  cart_info.read_geojson();

  if (!args.visual_file_name.empty()) {

    // Read visual variables (e.g., area and color) from CSV
    cart_info.read_csv();
  }

  // Project to equal area, if necessary
  // Write Output and EXIT if --output_equal_area_map is set to true
  cart_info.project_to_equal_area();

  // Store total number of GeoDivs to monitor progress
  size_t total_geo_divs = cart_info.n_geo_divs();

  // Write input map, with insets nicely placed
  if (args.plot_polygons) {
    cart_info.plot_input();
  }

  // Move insets to desired positions and EXIT
  if (args.output_shifted_insets) {
    cart_info.write_shifted_insets();
  }

  // Track progress of the cartogram generation
  ProgressTracker progress_tracker(static_cast<double>(total_geo_divs));

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

  for (InsetState &inset_state : cart_info.ref_to_inset_states()) {
    std::vector<std::pair<size_t, size_t>> changes =
      inset_state.densification_changes;
    double max_ratio = 0.0;
    double avg_ratio = 0.0;
    std::pair<size_t, size_t> worst_result = {0, 0};
    for (size_t i = 0; i < changes.size(); ++i) {
      double ratio = static_cast<double>(changes[i].second) / changes[i].first;
      std::cerr << "Iteration: " << i + 1 << ", "
                << "From: " << changes[i].first << ", "
                << "To: " << changes[i].second << ", "
                << "Ratio: " << ratio << std::endl;
      avg_ratio += ratio;
      if (ratio > max_ratio) {
        max_ratio = ratio;
        worst_result = changes[i];
      }
    }
    avg_ratio /= changes.size();
    std::cerr << "Maximum ratio: " << max_ratio << std::endl;
    std::cerr << "Average ratio: " << avg_ratio << std::endl;
    std::string csv_file_name = "densification_changes.csv";
    if (!std::filesystem::exists(csv_file_name)) {
      std::ofstream out_file_csv(csv_file_name);
      if (!out_file_csv) {
        std::cerr << "ERROR writing CSV: failed to open " << csv_file_name
                  << std::endl;
      } else {
        out_file_csv << "csv_file,from,to,max_ratio,avg_ratio\n";
      }
    }

    // Append to the CSV file
    std::ofstream out_file_csv(csv_file_name, std::ios_base::app);
    out_file_csv << inset_state.inset_name() << "," << worst_result.first
                 << "," << worst_result.second << ","
                 << static_cast<double>(worst_result.second) /
                      worst_result.first
                 << "," << avg_ratio << "\n";
    out_file_csv.close();
  }

  // Stop total time timer, and print time summary report
  // Export report with area errors to CSV if requested
  cart_info.print_time_report();
  if (!cart_info.converged())
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
