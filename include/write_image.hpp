#ifndef WRITE_IMAGE_H
#define WRITE_IMAGE_H

#include "colors.hpp"
#include "geo_div.hpp"
#include "inset_state.hpp"
#include <cairo/cairo-svg.h>

// TODO: IS THERE A CGAL WAY OF DETERMINING WHETHER THE LABEL'S BOUNDING
//       BOX IS COMPLETELY CONTAINED IN THE POLYGON?

// TODO: SHOULD THE CRITERION FOR PRINTING A LABEL BE THAT IT FITS INSIDE THE
//       POLYGON WITH HOLES? THAT CRITERION WOULD BE MORE RESTRICTIVE THAN
//       FITTING INSIDE THE EXTERIOR RING.

// ======================== Basic Plotting ========================

void write_polygon_points_on_surface(
  cairo_t *,
  Color,
  InsetState &inset_state);

void write_quadtree_rectangles_on_surface(
  cairo_t *cr,
  const std::vector<Bbox> &quadtree_bboxes,
  const Color &clr,
  const unsigned int ly);

void add_white_background(cairo_t *cr);

bool all_points_inside_exterior_ring(
  const std::vector<Point> &pts,
  const Polygon_with_holes &pwh);

void write_polygons_on_surface(
  cairo_t *cr,
  const bool fill_polygons,
  const bool colors,
  const InsetState &inset_state,
  const double line_width = 0.0,
  const Color clr = Color{0.0, 0.0, 0.0});

void write_grid_on_surface(cairo_t *cr, const InsetState &inset_state);

void write_cairo_polygons_to_svg(
  const std::string &,
  bool,
  bool,
  bool,
  bool,
  const std::unordered_map<Point, Vector> &,
  const InsetState &inset_state);

void write_labels_on_surface(cairo_t *cr, InsetState &inset_state);

void write_cells_on_surface(cairo_t *cr, InsetState &inset_state);

// ======================== Grid Cell ========================

double grid_cell_area(
  unsigned int x,
  unsigned int y,
  unsigned int cell_width,
  InsetState &inset_state);

std::pair<double, double> max_and_min_grid_cell_area(
  unsigned int cell_width,
  InsetState &inset_state);

std::pair<Point, Point> max_and_min_grid_cell_area_index(
  unsigned int cell_width,
  InsetState &inset_state);

double grid_cell_target_area(
  const unsigned int i,
  const unsigned int j,
  const double total_target_area,
  const double total_inset_area,
  InsetState &inset_state);

double grid_cell_target_area_per_km(
  const unsigned int i,
  const unsigned int j,
  const double total_target_area,
  const double total_inset_area,
  InsetState &inset_state);

double grid_cell_area_km(
  const InsetState &inset_state,
  const unsigned int i = 0,
  const unsigned int j = 0);

double compute_per_grid_cell(const InsetState &inset_state);

// ======================== Grid Heatmap ========================

void write_grid_heatmap_image(
  const std::string filename,
  const bool plot_equal_area_map,
  const bool crop_polygons,
  InsetState &inset_state);

std::vector<int> get_nice_numbers_for_bar(double max_target_area_per_km);

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

Bbox get_bbox_bar(
  const double bar_width,
  const double bar_height,
  InsetState &inset_state);

// ======================== Image Color ========================

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

std::vector<std::vector<Color>> grid_cell_colors(
  unsigned int cell_width,
  InsetState &inset_state);

// ======================== Legend Plotting ========================

void write_legend_on_surface(
  cairo_t *cr,
  bool equal_area_map,
  const InsetState &inset_state);

std::pair<unsigned int, unsigned int> get_km_legend_length(
  const InsetState &inset_state);

std::pair<unsigned int, unsigned int> get_visual_variable_legend_length(
  const InsetState &inset_state);

std::pair<std::string, std::string> get_legend_labels(
  unsigned int,
  unsigned int);

int get_nearest_nice_number_for_legend(int value);

Polygon transform_to_equal_area_projection_coor(
  Polygon edge_points,
  const InsetState &inset_state);

// Given area in the equal_area_projection projection coordinate system,
// returns the corresponding area in the square km^2
double equal_area_projection_area_to_earth_area(
  double equal_area_projection_area);

double earth_area_to_equal_area_projection_area(const double earth_area);

// ======================== Util ========================

double font_size(
  cairo_t *cr,
  const char *label,
  const Point label_pt,
  const GeoDiv &gd);

double get_font_size(
  cairo_t *cr,
  const char *label,
  const Point label_pt,
  const GeoDiv gd);

// Prints right-aligned text at specified location
void right_aligned_text(
  cairo_t *cr,
  double number,
  double x,
  double y,
  double font_size);

#endif  // WRITE_IMAGE_H
