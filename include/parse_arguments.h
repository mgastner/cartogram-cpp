#ifndef PARSE_ARGUMENTS_H_
#define PARSE_ARGUMENTS_H_

#include <iostream>
#include "argparse.hpp"

// Function to parse arguments and set variables in main()
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
  bool &simplification,
  bool &make_csv,
  bool &output_equal_area,
  bool &output_to_stdout,
  bool &plot_density,
  bool &plot_graticule,
  bool &plot_intersections,
  bool &plot_polygons);

#endif
