#ifndef CONST_H_
#define CONST_H_

#include <limits>
#include <numbers>

constexpr unsigned int default_long_grid_length = 512;
constexpr double dbl_epsilon = std::numeric_limits<double>::epsilon();
constexpr double dbl_inf = std::numeric_limits<double>::infinity();
constexpr double dbl_resolution = 1e-8;
constexpr unsigned int max_integrations = 100;
constexpr double max_permitted_area_error = 0.01;
constexpr double padding_unless_world = 1.5;
constexpr double pi = std::numbers::pi;
constexpr double earth_surface_area = 510.1e6;

// The resolution represents the number of rays to shoot through each cell
// (see inset_state/scanline_graph.cpp).
// default_resolution is used for filling grid cells with density.
// intersection_resolution is used to specify the number of scanlines shown
// in output files showing intersections.
constexpr unsigned int default_resolution = 16;
constexpr unsigned int intersections_resolution = 1;

// Points after simplification
constexpr unsigned int default_target_points_per_inset = 10000;
constexpr unsigned int min_points_per_ring = 10;

// Fraction of square side length by which squares on heatmap overlap
constexpr double sq_overlap = 0.2;

// Percent of total width/height (whichever is greater) of all insets that
// should be empty space
constexpr double inset_spacing_factor = 0.1;

// Percent of height/width of tallest/widest inset that divider should be
constexpr double divider_length = 0.8;

// Font size range for cairo labelling
constexpr double min_font_size = 6.0;
constexpr double max_font_size = 10.0;

// Threshold as a fraction of non-na and non-zero total area for a target
// area to be considered "too small"
constexpr double small_area_threshold_frac = 2e-5;

#endif
