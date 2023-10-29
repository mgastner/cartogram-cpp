#ifndef WRITE_IMAGE_H
#define WRITE_IMAGE_H

#include "colors.h"
#include "constants.h"
#include "inset_state.h"
#include <cairo/cairo-svg.h>

// TODO: IS THERE A CGAL WAY OF DETERMINING WHETHER THE LABEL'S BOUNDING
//       BOX IS COMPLETELY CONTAINED IN THE POLYGON?

// TODO: SHOULD THE CRITERION FOR PRINTING A LABEL BE THAT IT FITS INSIDE THE
//       POLYGON WITH HOLES? THAT CRITERION WOULD BE MORE RESTRICTIVE THAN
//       FITTING INSIDE THE EXTERIOR RING.
bool all_points_inside_exterior_ring(
  const std::vector<Point> &pts,
  const Polygon_with_holes &pwh);

double get_font_size(
  cairo_t *cr,
  const char *label,
  const Point label_pt,
  const GeoDiv gd);

// Functions deal with colors
Color interpolate_color(
  const double x,
  const double xmin,
  const double xmax,
  const Color ymin,
  const Color ymax);

// Diverging colour palette, mean accounted for
Color heatmap_color(
  const double dens,
  const double dens_min,
  const double dens_mean,
  const double dens_max);

// Sequential colour palette, mean not accounted for
Color grid_cell_color(
  const double area,
  const double max_area,
  const double min_area);

// Prints right-aligned text at specified location
void right_aligned_text(
  cairo_t *cr,
  double number,
  double x,
  double y,
  double font_size);

// Writes heatmap bar to cairo surface
void write_grid_heatmap_bar_to_cairo_surface(
  double min_value,
  double max_value,
  cairo_t *cr,
  Bbox bbox_bar,
  std::vector<std::pair<double, double>> major_ticks,
  std::vector<std::pair<double, double>> minor_ticks,
  const unsigned int ly);

// Given area in the albers projection coordinate system, returns the
// corresponding area in the square km^2
double albers_area_to_earth_area(const double albers_area);

std::vector<std::pair<double, double>> get_major_ticks(
  const double min_target_area_per_km,
  const double max_target_area_per_km,
  const double min_area_cell_point_area,
  const double max_area_cell_point_area,
  std::vector<int> nice_numbers);

std::vector<std::pair<double, double>> get_minor_ticks(
  int n_ticks_per_major,
  const double min_target_area_per_km,
  const double max_target_area_per_km,
  const double min_area_cell_point_area,
  const double max_area_cell_point_area,
  std::vector<int> nice_numbers);

std::pair<
  std::vector<std::pair<double, double>>,
  std::vector<std::pair<double, double>>>
get_ticks(
  const int n_ticks_per_major,
  const double min_target_area_per_km,
  const double max_target_area_per_km,
  const double min_area_cell_point_area,
  const double max_area_cell_point_area,
  std::vector<int> nice_numbers);

std::vector<int> get_nice_numbers_for_bar(const double max_target_area_per_km);

#endif
