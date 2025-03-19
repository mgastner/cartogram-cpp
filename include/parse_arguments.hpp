#ifndef PARSE_ARGUMENTS_HPP_
#define PARSE_ARGUMENTS_HPP_

#include "argparse.hpp"

struct Arguments {

  // Name of Geometry file (GeoJSON) and Visual Variables file (CSV)
  std::string geo_file_name;
  std::string visual_file_name;

  // Default number of grid cells along longer Cartesian coordinate axis
  unsigned int n_grid_rows_or_cols;

  // Target number of points to retain after simplification
  unsigned int target_points_per_inset;

  // Minimum number of integrations regardless of area error reached
  unsigned int min_integrations;

  // World maps need special projections
  bool world;

  // If `triangulation` is true, we apply a cartogram projection method based
  // on the triangulation of grid cells. It can eliminate intersections
  // that occur when the projected grid lines are strongly curved. Only
  // use this method if the tracer points are an FTReal2d data structure.
  bool triangulation;

  // Use Quadtree-Delaunay triangulation method
  bool qtdt_method;

// Should the polygons be simplified and densified?
  bool simplify;

  // If `rays` is true, we use the ray-shooting method to fill the grid cells.
  bool rays;

  // output_* writes said output and exits
  bool output_equal_area_map;
  bool output_shifted_insets;

  // export_* writes said outputs, and continues running
  bool redirect_exports_to_stdout;
  bool export_preprocessed;
  bool export_time_report;

  // Remove tiny polygons below threshold?
  // Criteria: If the proportion of the polygon area is smaller than
  // min_polygon_area * total area, then remove polygon
  bool remove_tiny_polygons;
  double min_polygon_area;

  // Other boolean values that are needed to parse the command line arguments
  bool make_csv;
  bool plot_density;
  bool plot_grid;
  bool plot_intersections;
  bool plot_polygons;
  bool plot_quadtree;
  bool skip_projection;

  // Column names in provided visual variables file (CSV)
  std::optional<std::string> id_col;
  std::optional<std::string> area_col;
  std::string inset_col;
  std::string color_col;
  std::string label_col;
};

// Function to parse arguments and set variables in main()
void parsed_arguments(
  int argc,
  const char *argv[],
  Arguments &args);

#endif // PARSE_ARGUMENTS_HPP_
