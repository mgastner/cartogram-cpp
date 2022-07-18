#ifndef PARSE_ARGUMENTS_H_
#define PARSE_ARGUMENTS_H_

#include "argparse.hpp"
#include <iostream>

// Function to parse arguments and set variables in main()
argparse::ArgumentParser parsed_arguments(
  int argc,
  const char *argv[],
  std::string &geo_file_name,
  std::string &visual_file_name,
  unsigned int &max_n_grid_rows_or_cols,
  unsigned int &target_points_per_inset,
  std::string &compare_geo_file_name,
  bool &insert_visual_variable,
  bool &world,
  bool &triangulation,
  bool &qtdt_method,
  bool &simplify,
  bool &make_csv,
  bool &output_equal_area,
  bool &output_to_stdout,
  bool &plot_density,
  bool &plot_graticule,
  bool &plot_intersections,
  bool &plot_polygons,
  bool &remove_tiny_polygons,
  double &minimum_polygon_area,
  bool &plot_quadtree);

#endif
