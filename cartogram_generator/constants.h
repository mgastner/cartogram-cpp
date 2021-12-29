#ifndef CONST_H_
#define CONST_H_

#include <limits>
#include <numbers>

constexpr unsigned int default_long_grid_side_length = 512;
constexpr double dbl_epsilon = std::numeric_limits<double>::epsilon();
constexpr double dbl_inf = std::numeric_limits<double>::infinity();
constexpr unsigned int max_integrations = 200;
constexpr double max_permitted_area_error = 0.005;
constexpr double padding_unless_world = 1.5;
constexpr double pi = std::numbers::pi;
constexpr double earth_radius = 6378.1370;

// Points after simplification
constexpr unsigned int target_points_per_inset = 10000;
constexpr unsigned int min_points_per_ring = 10;

// Fraction of square side length by which squares on heatmap overlap
constexpr double sq_overlap = 0.2;

// Percent of total width/height (whichever is greater) of all insets that
// should be empty space
constexpr double inset_spacing_factor = 0.1;

// Percent of height/width of tallest/widest inset that divider should be
constexpr double divider_length = 0.8;

#endif
