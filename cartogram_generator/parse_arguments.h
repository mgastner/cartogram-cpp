#ifndef PARSE_ARGUMENTS_H_
#define PARSE_ARGUMENTS_H_

#include <iostream>
#include "argparse.hpp"

// Function to parse arguments and set variables in main()
argparse::ArgumentParser parsed_arguments(const int argc,
                                         const char *argv[],
                                         std::string &geo_file_name,
                                         std::string &visual_file_name,
                                         unsigned int &long_grid_side_length,
                                         bool &world,
                                         bool &triangulation,
                                         bool &make_csv,
                                         bool &make_polygon_eps,
                                         bool &output_equal_area,
                                         bool &output_to_stdout,
                                         bool &plot_density,
                                         bool &plot_graticule);

#endif
