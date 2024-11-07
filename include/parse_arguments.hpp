#ifndef PARSE_ARGUMENTS_HPP_
#define PARSE_ARGUMENTS_HPP_

#include "argparse.hpp"

// Function to parse arguments and set variables in main()
argparse::ArgumentParser parsed_arguments(
  int argc,
  const char *argv[],
  std::string &geo_file_name,
  std::string &visual_file_name,
  unsigned int &max_n_grid_rows_or_cols,
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
  bool &rays,
  bool &output_preprocessed);

#endif // PARSE_ARGUMENTS_HPP_
