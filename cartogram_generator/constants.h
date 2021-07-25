#ifndef CONST_H_
#define CONST_H_

#include <limits>
#include <numbers>

constexpr unsigned int default_long_grid_side_length = 512;
constexpr double dbl_inf = std::numeric_limits<double>::infinity();
constexpr unsigned int max_integrations = 100;
constexpr double max_permitted_area_error = 0.01;
constexpr double padding_unless_world = 1.5;
constexpr double pi = std::numbers::pi;

// Fraction of square side length by which squares on heatmap overlap
constexpr double sq_overlap = 0.2;

#endif
