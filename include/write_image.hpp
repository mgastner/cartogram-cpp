#ifndef WRITE_IMAGE_H
#define WRITE_IMAGE_H

#include "colors.hpp"
#include "geo_div.hpp"
#include <cairo/cairo-svg.h>

// TODO: IS THERE A CGAL WAY OF DETERMINING WHETHER THE LABEL'S BOUNDING
//       BOX IS COMPLETELY CONTAINED IN THE POLYGON?

// TODO: SHOULD THE CRITERION FOR PRINTING A LABEL BE THAT IT FITS INSIDE THE
//       POLYGON WITH HOLES? THAT CRITERION WOULD BE MORE RESTRICTIVE THAN
//       FITTING INSIDE THE EXTERIOR RING.
bool all_points_inside_exterior_ring(
  std::vector<Point> &pts,
  Polygon_with_holes &pwh);

double get_font_size(cairo_t *cr, char *label, Point label_pt, GeoDiv gd);

// Functions deal with colors
Color interpolate_color(
  double x,
  double xmin,
  double xmax,
  Color ymin,
  Color ymax);

// Diverging colour palette, mean accounted for
Color heatmap_color(
  double dens,
  double dens_min,
  double dens_mean,
  double dens_max);

// Sequential colour palette, mean not accounted for
Color grid_cell_color(double area, double max_area, double min_area);

// Prints right-aligned text at specified location
void right_aligned_text(
  cairo_t *cr,
  double number,
  double x,
  double y,
  double font_size);

// Given area in the equal_area_projection projection coordinate system,
// returns the corresponding area in the square km^2
double equal_area_projection_area_to_earth_area(
  double equal_area_projection_area);
double earth_area_to_equal_area_projection_area(
  const double earth_area);

std::vector<std::pair<double, double>> get_major_ticks(
  double min_target_area_per_km,
  double max_target_area_per_km,
  double min_area_cell_point_area,
  double max_area_cell_point_area,
  std::vector<int> nice_numbers);

std::vector<std::pair<double, double>> get_minor_ticks(
  int n_ticks_per_major,
  double min_target_area_per_km,
  double max_target_area_per_km,
  double min_area_cell_point_area,
  double max_area_cell_point_area,
  std::vector<int> nice_numbers);

std::pair<
  std::vector<std::pair<double, double>>,
  std::vector<std::pair<double, double>>>
get_ticks(
  int n_ticks_per_major,
  double min_target_area_per_km,
  double max_target_area_per_km,
  double min_area_cell_point_area,
  double max_area_cell_point_area,
  std::vector<int> nice_numbers);

std::vector<int> get_nice_numbers_for_bar(double max_target_area_per_km);

#endif // WRITE_IMAGE_H
