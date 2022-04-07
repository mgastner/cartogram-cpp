#include "cartogram_info.h"
#include "colors.h"
#include "constants.h"
#include "inset_state.h"

#include <array>
#include <cairo/cairo-ps.h>
#include <cairo/cairo-svg.h>
#include <cairo/cairo.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void write_ps_header(const std::string filename, cairo_surface_t *surface)
{
  const std::string title = "%%Title: " + filename;
  cairo_ps_surface_dsc_comment(surface, title.c_str());
  cairo_ps_surface_dsc_comment(surface, "%%Creator: Michael T. Gastner et al.");
  cairo_ps_surface_dsc_comment(surface, "%%For: Humanity");
  cairo_ps_surface_dsc_comment(surface, "%%Copyright: License CC BY");
  cairo_ps_surface_dsc_comment(surface, "%%Magnification: 1.0000");
}

// TODO: IS THERE A CGAL WAY OF DETERMINING WHETHER THE LABEL'S BOUNDING
//       BOX IS COMPLETELY CONTAINED IN THE POLYGON?

// TODO: SHOULD THE CRITERION FOR PRINTING A LABEL BE THAT IT FITS INSIDE THE
//       POLYGON WITH HOLES? THAT CRITERION WOULD BE MORE RESTRICTIVE THAN
//       FITTING INSIDE THE EXTERIOR RING.
bool all_points_inside_exterior_ring(
  const std::vector<Point> pts,
  const Polygon_with_holes pwh)
{
  for (const auto &pt : pts) {
    if (pwh.outer_boundary().has_on_unbounded_side(pt)) {
      return false;
    }
  }
  return true;
}

double
font_size(cairo_t *cr, const char *label, const Point label_pt, const GeoDiv gd)
{
  for (double fsize = max_font_size; fsize >= min_font_size; fsize -= 0.5) {
    cairo_set_font_size(cr, fsize);
    cairo_text_extents_t extents;
    cairo_text_extents(cr, label, &extents);
    const Polygon_with_holes largest_pwh = gd.largest_polygon_with_holes();

    // Bounding box of the label
    const CGAL::Bbox_2 bb(
      label_pt.x() - 0.5 * extents.width,
      label_pt.y() - 0.5 * extents.height,
      label_pt.x() + 0.5 * extents.width,
      label_pt.y() + 0.5 * extents.height);

    // Vector of bounding-box edge points
    std::vector<Point> bb_edge_points;
    for (unsigned int i = 0; i <= 1; ++i) {
      for (unsigned int j = 0; j <= 5; ++j) {
        bb_edge_points.push_back(Point(
          (j * bb.xmin() + (5 - j) * bb.xmax()) / 5,
          (i * bb.ymin() + (1 - i) * bb.ymax())));
      }
    }
    if (all_points_inside_exterior_ring(bb_edge_points, largest_pwh)) {
      return fsize;
    }
  }
  return 0.0;
}

void InsetState::write_labels_to_cairo_surface(cairo_t *cr)
{
  for (const auto &gd : geo_divs()) {
    const std::string label = label_at(gd.id());
    const char *const label_char = label.c_str();
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_select_font_face(
      cr,
      "sans-serif",
      CAIRO_FONT_SLANT_NORMAL,
      CAIRO_FONT_WEIGHT_NORMAL);
    const Point label_pt = gd.point_on_surface_of_geodiv();

    // Get size of label
    const double fsize = font_size(cr, label_char, label_pt, gd);
    cairo_text_extents_t extents;

    // Draw label only if appropriate size is found
    if (fsize > 0.0) {
      cairo_set_font_size(cr, fsize);
      cairo_text_extents(cr, label_char, &extents);
      const double x = label_pt.x() - (extents.width / 2 + extents.x_bearing);
      const double y =
        ly_ - label_pt.y() - (extents.height / 2 + extents.y_bearing);
      cairo_move_to(cr, x, y);
      cairo_show_text(cr, label_char);
    }
  }
}

// Returns edge points of a graticule cell and treats
// them as a polygon.
Polygon InsetState::graticule_cell_edge_points(
  unsigned int x,
  unsigned int y,
  unsigned int cell_width = graticule_width,
  bool plot_equal_area_map = false)
{
  Polygon cell_edge_points;
  const boost::multi_array<XYPoint, 2> &proj =
    plot_equal_area_map ? original_proj_ : cum_proj_;

  // Horizontal lower edge points
  for (unsigned int i = x; i < x + cell_width; ++i) {
    double x_coor_trans = proj[i][y].x;
    double y_coor_trans = proj[i][y].y;
    cell_edge_points.push_back(Point(x_coor_trans, y_coor_trans));
  }

  // Vertical right edge points
  for (unsigned int i = y; i < y + cell_width; ++i) {
    double x_coor_trans = proj[x + cell_width][i].x;
    double y_coor_trans = proj[x + cell_width][i].y;
    cell_edge_points.push_back(Point(x_coor_trans, y_coor_trans));
  }

  // Horizontal upper edge points
  for (unsigned int i = x + cell_width; i > x; --i) {
    double x_coor_trans = proj[i][y + cell_width].x;
    double y_coor_trans = proj[i][y + cell_width].y;
    cell_edge_points.push_back(Point(x_coor_trans, y_coor_trans));
  }

  // Vertical left edge points
  for (unsigned int i = y + cell_width; i > y; --i) {
    double x_coor_trans = proj[x][i].x;
    double y_coor_trans = proj[x][i].y;
    cell_edge_points.push_back(Point(x_coor_trans, y_coor_trans));
  }

  // Complete the polygon by making first and last point same
  cell_edge_points.push_back(cell_edge_points[0]);

  return cell_edge_points;
}

// Returns graticule cell area based on edge points
double InsetState::graticule_cell_area(
  unsigned int x,
  unsigned int y,
  unsigned int cell_width)
{
  // Taking absolule to ensure we get the area irrespective of direction
  return abs(graticule_cell_edge_points(x, y, cell_width).area());
}

// Returns the largest and smallest graticule cell area to be used for
// graticule heatmap generation
std::pair<double, double>
InsetState::max_and_min_graticule_cell_area(unsigned int cell_width)
{
  // Initialize max and min area
  double max_area = -dbl_inf;
  double min_area = dbl_inf;

  // Iterate over graticule cells
  for (unsigned int i = 0; i < lx_ - cell_width; i += cell_width) {
    for (unsigned int j = 0; j < ly_ - cell_width; j += cell_width) {
      const double area = graticule_cell_area(i, j, cell_width);
      max_area = std::max(max_area, area);
      min_area = std::min(min_area, area);
    }
  }
  return std::make_pair(max_area, min_area);
}

std::pair<Point, Point>
InsetState::max_and_min_graticule_cell_area_index(unsigned int cell_width)
{

  // Initialize max and min area
  double max_area = -dbl_inf;
  double min_area = dbl_inf;
  unsigned int max_i = 0;
  unsigned int max_j = 0;
  unsigned int min_i = 0;
  unsigned int min_j = 0;

  // Iterate over graticule cells
  for (unsigned int i = 0; i < lx_ - cell_width; i += cell_width) {
    for (unsigned int j = 0; j < ly_ - cell_width; j += cell_width) {
      const double area = graticule_cell_area(i, j, cell_width);
      if (area > max_area) {
        max_area = area;
        max_i = i;
        max_j = j;
      }
      if (area < min_area) {
        min_area = area;
        min_i = i;
        min_j = j;
      }
    }
  }
  return std::make_pair(Point(max_i, max_j), Point(min_i, min_j));
}

// Functions deal with colors
Color interpolate_color(
  const double x,
  const double xmin,
  const double xmax,
  const Color ymin,
  const Color ymax)
{
  Color interpolated_color;
  // Interpolate color for red, green and blue value
  for (char c : {'r', 'g', 'b'}) {
    interpolated_color(c) =
      ((x - xmin) * ymax(c) + (xmax - x) * ymin(c)) / (xmax - xmin);
  }
  return interpolated_color;
}

// Diverging colour palette, mean accounted for
Color heatmap_color(
  const double dens,
  const double dens_min,
  const double dens_mean,
  const double dens_max)
{

  // Assign possible categories for red, green, blue
  const std::vector<Color> colors = {

    // Red to blue
    Color("#67001f"),
    Color("#b2182b"),
    Color("#d6604d"),
    Color("#f4a582"),
    Color("#fddbc7"),
    Color("#f7f7f7"),
    Color("#d1e5f0"),
    Color("#92c5de"),
    Color("#4393c3"),
    Color("#2166ac"),
    Color("#053061")

    // // Turqoise to brown
    // Color("#543005"),
    // Color("#8c510a"),
    // Color("#bf812d"),
    // Color("#dfc27d"),
    // Color("#f6e8c3"),
    // Color("#f5f5f5"),
    // Color("#c7eae5"),
    // Color("#80cdc1"),
    // Color("#35978f"),
    // Color("#01665e"),
    // Color("#003c30")

    // // Original
    // Color(0.33, 0.19, 0.02),
    // Color(0.55, 0.32, 0.04),
    // Color(0.75, 0.51, 0.18),
    // Color(0.87, 0.76, 0.49),
    // Color(0.96, 0.91, 0.76),
    // Color(0.99, 0.96, 0.89),
    // Color(0.78, 0.92, 0.90),
    // Color(0.50, 0.80, 0.76),
    // Color(0.21, 0.59, 0.56),
    // Color(0.00, 0.40, 0.37),
    // Color(0.00, 0.24, 0.19)
  };
  int n_categories = colors.size();
  double xmin, xmax;
  int color_category;

  // If no discernible difference between dens and miniimum density, set
  // lowest
  if (std::fabs(dens - dens_min) <= dbl_resolution) {
    return colors[n_categories - 1];
  }

  // Choose color category
  if (dens >= dens_max) {
    return colors[0];
  } else if (dens > dens_mean) {
    color_category = 5 * (dens_max - dens) / (dens_max - dens_mean);
    xmax = dens_max - 0.2 * color_category * (dens_max - dens_mean);
    xmin = xmax - 0.2 * (dens_max - dens_mean);

    // Assign color category 0 if dens_max and dens are very close
    color_category = std::max(color_category, 0);
  } else if (dens > dens_min) {
    color_category = 5 * (dens_mean - dens) / (dens_mean - dens_min) + 5;
    xmax = dens_mean - 0.2 * (color_category - 5) * (dens_mean - dens_min);
    xmin = xmax - 0.2 * (dens_mean - dens_min);

    // Assign color category 9 if dens_min and dens are very close
    color_category = std::min(color_category, n_categories - 2);
  } else {
    return colors[n_categories - 1];
  }
  return interpolate_color(
    dens,
    xmin,
    xmax,
    colors[color_category + 1],
    colors[color_category]);
}

// Sequential colour palette, mean not accounted for
Color graticule_cell_color(
  const double area,
  const double max_area,
  const double min_area)
{
  // Assign possible categories for red, green, blue
  const std::vector<Color> colors = {

    // White to purple
    Color("#fcfbfd"),
    Color("#efedf5"),
    Color("#dadaeb"),
    Color("#bcbddc"),
    Color("#9e9ac8"),
    Color("#807dba"),
    Color("#6a51a3"),
    Color("#54278f"),
    Color("#3f007d")

    // White to green
    // Color("#f7fcf5"),
    // Color("#e5f5e0"),
    // Color("#c7e9c0"),
    // Color("#a1d99b"),
    // Color("#74c476"),
    // Color("#41ab5d"),
    // Color("#238b45"),
    // Color("#006d2c"),
    // Color("#00441b")

    // // White to red
    // Color(1.000, 0.961, 0.941),
    // Color(0.996, 0.878, 0.824),
    // Color(0.988, 0.733, 0.631),
    // Color(0.988, 0.572, 0.447),
    // Color(0.984, 0.416, 0.290),
    // Color(0.937, 0.231, 0.173),
    // Color(0.796, 0.094, 0.114),
    // Color(0.647, 0.058, 0.082),
    // Color(0.404, 0.000, 0.050)
  };
  int n_categories = colors.size();

  // Normalize area to [0,1] and make it logarithmic
  double ratio = (log(area) - log(min_area)) / (log(max_area) - log(min_area));

  // Determine color category
  double category = fmax(0, ratio * (n_categories - 1) - 10e-6);
  int xmin = floor(category);
  int xmax = ceil(category);

  if (area == max_area) {
    return colors[n_categories - 1];
  } else if (area == min_area) {
    return colors[0];
  } else {
    Color interpolated_color;
    for (char c : {'r', 'g', 'b'}) {
      interpolated_color(c) =
        colors[xmin](c) +
        (colors[xmax](c) - colors[xmin](c)) * (category - xmin);
    }
    return interpolated_color;
  }
}

std::vector<std::vector<Color>>
InsetState::graticule_cell_colors(unsigned int cell_width)
{

  // Initialize max and min area
  double max_area, min_area;
  std::tie(max_area, min_area) = max_and_min_graticule_cell_area(cell_width);

  // Initialize colors
  std::vector<std::vector<Color>> colors(
    lx_ - cell_width,
    std::vector<Color>(ly_ - cell_width));

  // Iterate over graticule cells
  for (unsigned int i = 0; i < lx_ - cell_width; i += cell_width) {
    for (unsigned int j = 0; j < ly_ - cell_width; j += cell_width) {
      const double area = graticule_cell_area(i, j, cell_width);
      colors[i][j] = graticule_cell_color(area, max_area, min_area);
    }
  }
  return colors;
}

void write_graticule_heatmap_bar_to_cairo_surface(
  double min_value,
  double max_value,
  cairo_t *cr,
  Bbox bbox_bar,
  std::vector<std::pair<double, double>> major_ticks,
  std::vector<std::pair<double, double>> minor_ticks,
  const unsigned int ly)
{
  const int n_gradient_bars = 500;

  // get bar coordinates
  const double xmin_bar = bbox_bar.xmin();
  const double xmax_bar = bbox_bar.xmax();
  const double ymin_bar = bbox_bar.ymin();
  const double ymax_bar = bbox_bar.ymax();

  const double bar_width = xmax_bar - xmin_bar;

  // calculate individual bar gradient segment property
  const double gradient_segment_height =
    (ymax_bar - ymin_bar) / n_gradient_bars;
  const double gradient_segment_value =
    (max_value - min_value) / n_gradient_bars;

  // Draw the outer bar lines
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_set_line_width(cr, 1.0);
  cairo_move_to(cr, xmin_bar, ly - ymin_bar);
  cairo_line_to(cr, xmax_bar, ly - ymin_bar);
  cairo_line_to(cr, xmax_bar, ly - ymax_bar);
  cairo_line_to(cr, xmin_bar, ly - ymax_bar);
  cairo_line_to(cr, xmin_bar, ly - ymin_bar);
  cairo_stroke(cr);

  // Draw the gradient segment rectangles
  double value_at_gradient_segment = min_value;

  for (double y = ymin_bar; y <= ymax_bar; y += gradient_segment_height) {
    Color color = graticule_cell_color(
      exp(value_at_gradient_segment),
      exp(max_value),
      exp(min_value));
    cairo_set_source_rgb(cr, color.r, color.g, color.b);
    cairo_rectangle(cr, xmin_bar, ly - y, bar_width, gradient_segment_height);
    cairo_fill(cr);
    value_at_gradient_segment += gradient_segment_value;
  }

  unsigned int font_size = 8;

  // Set font properties
  cairo_select_font_face(
    cr,
    "Sans",
    CAIRO_FONT_SLANT_NORMAL,
    CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size(cr, font_size);

  // Draw the ticks and nice_numbers
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_set_line_width(cr, .9);
  for (auto tick : major_ticks) {
    double area = tick.first;
    int NiceNumber = tick.second;
    if (area > min_value and area < max_value) {
      double y =
        ((log(area) - log(min_value)) / (log(max_value) - log(min_value))) *
          (ymax_bar - ymin_bar) +
        ymin_bar;
      cairo_move_to(cr, xmax_bar - 6, ly - y);
      cairo_line_to(cr, xmax_bar, ly - y);
      cairo_move_to(cr, xmax_bar + 1, ly - y + font_size / 2.0);
      cairo_show_text(cr, std::to_string(NiceNumber).c_str());
      cairo_stroke(cr);
    }
  }
  cairo_set_line_width(cr, .7);
  for (auto ticks : minor_ticks) {
    double area = ticks.first;
    if (area > min_value and area < max_value) {
      double y =
        ((log(area) - log(min_value)) / (log(max_value) - log(min_value))) *
          (ymax_bar - ymin_bar) +
        ymin_bar;

      cairo_move_to(cr, xmax_bar - 3, ly - y);
      cairo_line_to(cr, xmax_bar, ly - y);
      cairo_stroke(cr);
    }
  }

  // TODO: Use cairo_text_extents_t for precise placement of text
  cairo_set_line_width(cr, 1.0);
  std::string bar_text_top = "Cases per";
  std::string bar_text_bottom = "km²";
  cairo_move_to(
    cr,
    (xmin_bar + xmax_bar) / 2 - 16,
    ly - ymax_bar - (font_size * 2.0));
  cairo_show_text(cr, bar_text_top.c_str());
  cairo_move_to(cr, (xmin_bar + xmax_bar) / 2 - 5.5, ly - ymax_bar - font_size);
  cairo_show_text(cr, bar_text_bottom.c_str());
}

void InsetState::trim_graticule_heatmap(cairo_t *cr, double padding)
{
  // Canvas dimension
  Bbox is_bbox = bbox();
  double xmin = is_bbox.xmin() - padding;
  double xmax = is_bbox.xmax() + padding;
  double ymin = is_bbox.ymin() - padding;
  double ymax = is_bbox.ymax() + padding;

  // Color white outside is_bbox
  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_rectangle(cr, 0, 0, xmin, ymax);
  cairo_rectangle(cr, xmax, 0, lx_ - xmax, ymax);
  cairo_rectangle(cr, 0, ymax, lx_, ly_ - ymax);
  cairo_rectangle(cr, 0, 0, lx_, ymin);
  cairo_fill(cr);
}

// Writes graticule cells and colors them if required
void InsetState::write_graticules_to_cairo_surface(cairo_t *cr)
{

  // Set line width of graticule lines
  cairo_set_line_width(cr, 5e-4 * std::min(lx_, ly_));
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);

  // Iterate over graticule cells
  for (unsigned int i = 0; i < lx_ - graticule_width; i += graticule_width) {
    for (unsigned int j = 0; j < ly_ - graticule_width; j += graticule_width) {
      // Draw graticule cell by connecting edge points
      const Polygon cell_edge_points = graticule_cell_edge_points(i, j);
      cairo_move_to(cr, cell_edge_points[0].x(), ly_ - cell_edge_points[0].y());
      for (unsigned int k = 1; k < cell_edge_points.size(); ++k) {
        cairo_line_to(
          cr,
          cell_edge_points[k].x(),
          ly_ - cell_edge_points[k].y());
      }
      cairo_stroke(cr);
    }
  }
}

void InsetState::write_graticule_colors_to_cairo_surface(
  cairo_t *cr,
  bool plot_equal_area_map,
  bool crop)
{
  unsigned int cell_width = 1;

  // Get colors
  const auto colors = graticule_cell_colors(cell_width);

  // Set line width of graticule lines
  cairo_set_line_width(cr, 5e-6 * std::min(lx_, ly_));

  // Print the color of all graticule cells and exit
  if (!crop) {

    // Iterate over graticule cells
    for (unsigned int i = 0; i < lx_ - cell_width; i += cell_width) {
      for (unsigned int j = 0; j < ly_ - cell_width; j += cell_width) {
        // Set color of the border of the graticule polygon
        cairo_set_source_rgb(
          cr,
          colors[i][j].r,
          colors[i][j].g,
          colors[i][j].b);

        // Draw graticule cell by connecting edge points
        const Polygon cell_edge_points =
          graticule_cell_edge_points(i, j, cell_width, plot_equal_area_map);
        cairo_move_to(
          cr,
          cell_edge_points[0].x(),
          ly_ - cell_edge_points[0].y());
        for (unsigned int k = 1; k < cell_edge_points.size(); ++k) {
          cairo_line_to(
            cr,
            cell_edge_points[k].x(),
            ly_ - cell_edge_points[k].y());
        }

        // Fill the graticule polygon with color
        cairo_fill_preserve(cr);
        cairo_stroke(cr);
      }
    }

    // Don't plot uncropped version
    return;
  }

  // Draw cartogram polygons or equal area map polygons
  const std::vector<GeoDiv> &geo_divs =
    plot_equal_area_map ? geo_divs_original_ : geo_divs_;

  const std::vector<GeoDiv> &geo_divs_original = geo_divs_original_;
  std::map<std::string, Bbox> original_bboxes;

  for (const auto &gd : geo_divs_original) {
    original_bboxes[gd.id()] = gd.bbox();
  }

  // Clip to shape
  for (const auto &gd : geo_divs) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto ext_ring = pwh.outer_boundary();
      cairo_new_path(cr);
      cairo_move_to(cr, ext_ring[0].x(), ly_ - ext_ring[0].y());

      // Plot each point in exterior ring
      for (unsigned int i = 1; i < ext_ring.size(); ++i) {
        cairo_line_to(cr, ext_ring[i].x(), ly_ - ext_ring[i].y());
      }

      // Close entire path
      cairo_close_path(cr);
      cairo_clip(cr);

      Bbox gd_bbox = original_bboxes.at(gd.id());

      // Iterate over graticule cells
      for (unsigned int i = gd_bbox.xmin() - 1; i < gd_bbox.xmax() + 1;
           i += cell_width) {
        for (unsigned int j = gd_bbox.ymin() - 1; j < gd_bbox.ymax() + 1;
             j += cell_width) {

          // Set color of the border of the graticule polygon
          cairo_set_source_rgb(
            cr,
            colors[i][j].r,
            colors[i][j].g,
            colors[i][j].b);

          // Draw graticule cells by connecting edge points
          const auto cell_edge_points =
            graticule_cell_edge_points(i, j, cell_width, plot_equal_area_map);

          // Move to first graticule edge
          cairo_move_to(
            cr,
            cell_edge_points[0].x(),
            ly_ - cell_edge_points[0].y());

          // Draw remaining graticule cells
          for (unsigned int k = 1; k < cell_edge_points.size(); ++k) {
            cairo_line_to(
              cr,
              cell_edge_points[k].x(),
              ly_ - cell_edge_points[k].y());
          }

          // Fill the graticule polygon with color
          cairo_fill_preserve(cr);
          cairo_stroke(cr);
        }
      }

      // Remove GeoDiv clip
      cairo_reset_clip(cr);
    }
  }
  return;
}

void InsetState::write_polygons_to_cairo_surface(
  cairo_t *cr,
  const bool fill_polygons,
  const bool colors,
  const bool plot_equal_area_map)
{
  cairo_set_line_width(cr, 1e-3 * std::min(lx_, ly_));

  // Draw cartogram polygons or equal area map polygons
  const std::vector<GeoDiv> &geo_divs =
    plot_equal_area_map ? geo_divs_original_ : geo_divs_;

  // Draw the shapes
  for (const auto &gd : geo_divs) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const Polygon ext_ring = pwh.outer_boundary();
      cairo_move_to(cr, ext_ring[0].x(), ly_ - ext_ring[0].y());

      // Plot each point in exterior ring
      for (unsigned int i = 1; i < ext_ring.size(); ++i) {
        cairo_line_to(cr, ext_ring[i].x(), ly_ - ext_ring[i].y());
      }

      // Close the exterior ring
      cairo_close_path(cr);

      // Draw holes
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        cairo_move_to(cr, (*h)[0].x(), ly_ - (*h)[0].y());
        const size_t hsize = (*h).size();
        for (unsigned int i = 1; i <= hsize; ++i) {
          cairo_line_to(cr, (*h)[i % hsize].x(), ly_ - (*h)[i % hsize].y());
        }
      }
      if (colors || fill_polygons) {
        if (is_input_target_area_missing(gd.id())) {
          // Fill path with dark gray
          cairo_set_source_rgb(cr, 0.9375, 0.9375, 0.9375);
        } else if (colors) {
          // Get color
          const Color col = color_at(gd.id());

          // Fill path
          cairo_set_source_rgb(cr, col.r, col.g, col.b);
        } else if (fill_polygons) {
          // Fill path with default color
          cairo_set_source_rgb(cr, 0.96, 0.92, 0.70);
        }
        cairo_fill_preserve(cr);
      }

      cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
      cairo_stroke(cr);
    }
  }
}

// Outputs a SVG/PS file with polygons, labels, and graticules (if required)
void InsetState::write_map_image(
  const std::string filename,
  const bool fill_polygons,
  const bool plot_graticule,
  const bool image_format_ps)
{

  // Check whether the map has all GeoDivs colored
  const bool colors = (colors_size() == n_geo_divs());
  cairo_surface_t *surface;

  // Create a cairo surface
  image_format_ps
    ? surface = cairo_ps_surface_create(filename.c_str(), lx_, ly_)
    : surface = cairo_svg_surface_create(filename.c_str(), lx_, ly_);
  cairo_t *cr = cairo_create(surface);

  // Write header
  if (image_format_ps) {
    write_ps_header(filename, surface);
  }

  // Draw polygons with color
  write_polygons_to_cairo_surface(cr, fill_polygons, colors, false);

  // Place labels
  write_labels_to_cairo_surface(cr);

  // Draw graticule without color
  if (plot_graticule) {
    write_graticules_to_cairo_surface(cr);
  }
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

Polygon InsetState::transform_to_albers_coor(Polygon edge_points)
{
  Transformation scale(CGAL::SCALING, latt_const_);

  Polygon cell_edge_points_albers = transform(scale, edge_points);

  return cell_edge_points_albers;
}

double albers_area_to_earth_area(const double albers_area)
{
  return (albers_area * earth_surface_area) / (4 * pi);
}

double
InsetState::graticule_cell_area_km(const unsigned int i, const unsigned int j)
{
  Polygon cell_edge_points = graticule_cell_edge_points(i, j, 1, true);
  const Polygon cell_edge_points_albers =
    transform_to_albers_coor(cell_edge_points);
  const double cell_area = cell_edge_points_albers.area();

  const double cell_area_km = albers_area_to_earth_area(cell_area);

  return cell_area_km;
}

double InsetState::graticule_cell_target_area(
  const unsigned int i,
  const unsigned int j,
  const double total_target_area,
  const double total_inset_area)
{
  const Polygon cell_edge_points = graticule_cell_edge_points(i, j, 1, false);
  const double cell_area = cell_edge_points.area();

  const double cell_target_area =
    (cell_area * total_target_area) / total_inset_area;

  return cell_target_area;
}

double InsetState::graticule_cell_target_area_per_km(
  const unsigned int i,
  const unsigned int j,
  const double total_target_area,
  const double total_inset_area)
{
  const double cell_target_area =
    graticule_cell_target_area(i, j, total_target_area, total_inset_area);
  const double cell_area_km = graticule_cell_area_km(i, j);

  const double cell_target_area_per_km = cell_target_area / cell_area_km;

  return cell_target_area_per_km;
}

std::vector<std::pair<double, double>> get_major_ticks(
  const double min_target_area_per_km,
  const double max_target_area_per_km,
  const double min_area_cell_point_area,
  const double max_area_cell_point_area,
  std::vector<int> nice_numbers)

{
  std::vector<std::pair<double, double>> ticks;
  for (auto niceNumber : nice_numbers) {
    double NiceNumberRatio = (niceNumber - min_target_area_per_km) /
                             (max_target_area_per_km - min_target_area_per_km);
    double area =
      min_area_cell_point_area +
      NiceNumberRatio * (max_area_cell_point_area - min_area_cell_point_area);
    ticks.push_back(std::make_pair(area, niceNumber));
  }
  return ticks;
}

std::vector<std::pair<double, double>> get_minor_ticks(
  int n_ticks_per_major,
  const double min_target_area_per_km,
  const double max_target_area_per_km,
  const double min_area_cell_point_area,
  const double max_area_cell_point_area,
  std::vector<int> nice_numbers)

{
  n_ticks_per_major += 2;
  std::vector<std::pair<double, double>> minor_ticks;
  const int n_major_ticks = nice_numbers.size();
  for (int i = 0; i < n_major_ticks - 1; i++) {
    double first_major_tick = (double)nice_numbers[i];
    double second_major_tick = (double)nice_numbers[i + 1];
    double minor_tick_ratio =
      (second_major_tick - first_major_tick) / (n_ticks_per_major - 1);
    for (int j = 1; j < n_ticks_per_major - 1; j++) {
      double minor_tick = first_major_tick + minor_tick_ratio * j;
      double NiceNumberRatio =
        (minor_tick - min_target_area_per_km) /
        (max_target_area_per_km - min_target_area_per_km);
      double area =
        min_area_cell_point_area +
        NiceNumberRatio * (max_area_cell_point_area - min_area_cell_point_area);
      minor_ticks.push_back(std::make_pair(area, minor_tick));
    }
  }
  return minor_ticks;
}

std::pair<
  std::vector<std::pair<double, double>>,
  std::vector<std::pair<double, double>>>
get_ticks(
  const int n_ticks_per_major,
  const double min_target_area_per_km,
  const double max_target_area_per_km,
  const double min_area_cell_point_area,
  const double max_area_cell_point_area,
  std::vector<int> nice_numbers)
{
  std::vector<std::pair<double, double>> major_ticks = get_major_ticks(
    min_target_area_per_km,
    max_target_area_per_km,
    min_area_cell_point_area,
    max_area_cell_point_area,
    nice_numbers);
  std::vector<std::pair<double, double>> minor_ticks = get_minor_ticks(
    n_ticks_per_major,
    min_target_area_per_km,
    max_target_area_per_km,
    min_area_cell_point_area,
    max_area_cell_point_area,
    nice_numbers);
  std::pair<
    std::vector<std::pair<double, double>>,
    std::vector<std::pair<double, double>>>
    ticks(major_ticks, minor_ticks);
  return ticks;
}

Bbox InsetState::get_bbox_bar(const double bar_width, const double bar_height)
{
  const Bbox is_bbox = bbox();

  // Position the bar 25 pixels to the right of the is_bbox
  const double xmin_bar = is_bbox.xmax() + 35;
  const double xmax_bar = xmin_bar + bar_width;

  // Position the bar at the middle of the is_bbox y coordinates
  double ymid_bar = (is_bbox.ymax() + is_bbox.ymin()) / 2 - 25;
  double ymin_bar = ymid_bar - bar_height / 2;
  double ymax_bar = ymid_bar + bar_height / 2;

  const Bbox bbox_bar(xmin_bar, ymin_bar, xmax_bar, ymax_bar);

  return bbox_bar;
}

std::vector<int> get_nice_numbers_for_bar(const double max_target_area_per_km)
{
  std::vector<int> nice_numbers;
  int NiceNumber = 1;
  nice_numbers.push_back(NiceNumber);
  while (NiceNumber < max_target_area_per_km) {
    NiceNumber = NiceNumber * 10;
    nice_numbers.push_back(NiceNumber);
  }
  return nice_numbers;
}

// Outputs a SVG/PS file of graticule heatmap
void InsetState::write_graticule_heatmap_image(
  const std::string filename,
  const bool plot_equal_area_map,
  const bool image_format_ps,
  const bool crop)
{
  // Whether to draw bar on the cairo surface
  const bool draw_bar = true;

  // Create a cairo surface
  cairo_surface_t *surface;
  image_format_ps
    ? surface = cairo_ps_surface_create(filename.c_str(), lx_, ly_)
    : surface = cairo_svg_surface_create(filename.c_str(), lx_, ly_);
  cairo_t *cr = cairo_create(surface);

  // Get inset areas
  const double total_ta = total_target_area();
  const double total_ia = total_inset_area();

  const Bbox bbox_bar = get_bbox_bar(15, 150);

  // Get the max and min graticule cell area points
  Point max_area_cell_point, min_area_cell_point;

  std::tie(max_area_cell_point, min_area_cell_point) =
    max_and_min_graticule_cell_area_index(1);

  const double max_area_cell_point_area =
    graticule_cell_area(max_area_cell_point.x(), max_area_cell_point.y(), 1);
  const double min_area_cell_point_area =
    graticule_cell_area(min_area_cell_point.x(), min_area_cell_point.y(), 1);

  // Get the max and min graticule cell target area per km
  double max_target_area_per_km = graticule_cell_target_area_per_km(
    max_area_cell_point.x(),
    max_area_cell_point.y(),
    total_ta,
    total_ia);
  double min_target_area_per_km = graticule_cell_target_area_per_km(
    min_area_cell_point.x(),
    min_area_cell_point.y(),
    total_ta,
    total_ia);
  std::cerr << std::endl;
  std::cerr << "Max target area per km: " << max_target_area_per_km
            << std::endl;
  std::cerr << "Min target area per km: " << min_target_area_per_km
            << std::endl;

  std::vector<int> nice_numbers =
    get_nice_numbers_for_bar(max_target_area_per_km);

  std::vector<std::pair<double, double>> major_ticks, minor_ticks;

  std::tie(major_ticks, minor_ticks) = get_ticks(
    10,
    min_target_area_per_km,
    max_target_area_per_km,
    min_area_cell_point_area,
    max_area_cell_point_area,
    nice_numbers);

  // Write header
  if (image_format_ps) {
    write_ps_header(filename, surface);
  }

  // Draw colors
  write_graticule_colors_to_cairo_surface(cr, plot_equal_area_map, crop);

  // Draw polygons without color
  write_polygons_to_cairo_surface(cr, false, false, plot_equal_area_map);

  trim_graticule_heatmap(cr, 20);

  write_graticule_heatmap_bar_to_cairo_surface(
    min_area_cell_point_area,
    max_area_cell_point_area,
    cr,
    bbox_bar,
    major_ticks,
    minor_ticks,
    ly_);

  if (draw_bar) {
    std::string bar_filename =
      "bar_" + filename.substr(0, filename.size() - 2) + "svg";
    // std::string bar_filename = "bar_" + filename;

    // Create a cairo bar_surface
    cairo_surface_t *bar_surface =
      cairo_svg_surface_create(bar_filename.c_str(), 160, 400);
    // cairo_surface_t *bar_surface;
    // image_format_ps
    //   ? bar_surface = cairo_ps_surface_create(bar_filename.c_str(), 160, 400)
    //   : bar_surface = cairo_svg_surface_create(bar_filename.c_str(), 160,
    //   400);
    cairo_t *bar_cr = cairo_create(bar_surface);

    // Write header
    // if (image_format_ps) {
    //   write_ps_header(bar_filename, bar_surface);
    // }

    // Write bar
    write_graticule_heatmap_bar_to_cairo_surface(
      min_area_cell_point_area,
      max_area_cell_point_area,
      bar_cr,
      Bbox(75.0, 200.0, 95.0, 350.0),
      major_ticks,
      minor_ticks,
      ly_);

    cairo_show_page(bar_cr);
    cairo_surface_destroy(bar_surface);
    cairo_destroy(bar_cr);
  }
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

// Function to show the density bar on the cairo surface
void write_density_bar_to_cairo_surface(
  const double min_value,
  const double mean_value,
  const double max_value,
  cairo_t *cr,
  Bbox bbox_bar,
  const unsigned int ly)
{
  const int n_gradient_bars = 500;

  // get bar coordinates
  const double xmin_bar = bbox_bar.xmin();
  const double xmax_bar = bbox_bar.xmax();
  const double ymin_bar = bbox_bar.ymin();
  const double ymax_bar = bbox_bar.ymax();

  const double bar_width = xmax_bar - xmin_bar;

  // position of mean line along bar
  const double ymean_bar = ((ymax_bar - ymin_bar) / (max_value - min_value)) *
                             (mean_value - min_value) +
                           ymin_bar;

  // calculate individual bar gradient segment property
  const double gradient_segment_height =
    (ymax_bar - ymin_bar) / n_gradient_bars;
  const double gradient_segment_value =
    abs(max_value - min_value) / n_gradient_bars;

  // Draw the outer bar lines
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_set_line_width(cr, 1.0);
  cairo_move_to(cr, xmin_bar, ly - ymin_bar);
  cairo_line_to(cr, xmax_bar, ly - ymin_bar);
  cairo_line_to(cr, xmax_bar, ly - ymax_bar);
  cairo_line_to(cr, xmin_bar, ly - ymax_bar);
  cairo_line_to(cr, xmin_bar, ly - ymin_bar);
  cairo_stroke(cr);

  // Draw the gradient segment rectangles
  double value_at_gradient_segment = min_value;

  for (double y = ymin_bar; y <= ymax_bar; y += gradient_segment_height) {
    Color color = heatmap_color(
      value_at_gradient_segment,
      min_value,
      mean_value,
      max_value);
    cairo_set_source_rgb(cr, color.r, color.g, color.b);
    cairo_rectangle(cr, xmin_bar, ly - y, bar_width, gradient_segment_height);
    cairo_fill(cr);
    value_at_gradient_segment += gradient_segment_value;
  }

  // Draw the mean line
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_set_line_width(cr, 0.9);
  cairo_move_to(cr, xmin_bar + bar_width / 2, ly - ymean_bar);
  cairo_line_to(cr, xmax_bar, ly - ymean_bar);
  cairo_stroke(cr);

  unsigned int font_size = 8;

  // Set font properties
  cairo_select_font_face(
    cr,
    "Sans",
    CAIRO_FONT_SLANT_NORMAL,
    CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size(cr, font_size);

  // Write "Density" on top of bar
  cairo_move_to(
    cr,
    xmin_bar - bar_width / 2 - 1,
    ly - ymax_bar - (font_size * 3.0));
  cairo_show_text(cr, "Residual");
  cairo_move_to(
    cr,
    xmin_bar - bar_width / 2 - 1,
    ly - ymax_bar - (font_size * 2.0));
  cairo_show_text(cr, "Density (km²)");

  std::string temp = std::to_string(max_value);
  // const double min_value, const double mean_value, const double max_value;
  // Write "High" top right of bar
  cairo_move_to(cr, xmax_bar + 5, ly - ymax_bar);
  cairo_show_text(cr, temp.substr(0, temp.size() - 4).c_str());

  // cairo_show_text(cr, "High");

  // Write "Low" bottom right of bar
  temp = std::to_string(min_value);
  cairo_move_to(cr, xmax_bar + 5, ly - ymin_bar + (font_size / 2.0));
  cairo_show_text(cr, temp.substr(0, temp.size() - 4).c_str());
  // cairo_show_text(cr, "Low");

  // Write "Mean" beside ymean_bar
  temp = std::to_string(mean_value);
  cairo_move_to(cr, xmax_bar + 5, ly - ymean_bar + (font_size / 4.0));
  cairo_show_text(cr, temp.substr(0, temp.size() - 4).c_str());

  font_size /= 2;

  // double max_absolute_val = std::max(std::abs(max_value),
  // std::abs(min_value)); Draw remaining ticks
  long long magnitude = std::pow(10, floor(std::log10(max_value)));

  if (max_value < magnitude * 2 && std::abs(min_value) < magnitude * 2) {
    magnitude /= 4;
  } else if (max_value < magnitude * 5 && std::abs(min_value) < magnitude * 5) {
    magnitude /= 2;
  }
  if (magnitude >= 0.5) {
    cairo_set_line_width(cr, 0.7);
    double bar_ratio = (ymax_bar - ymin_bar) / (max_value - min_value);

    // Positive ticks
    for (unsigned int i = magnitude; i < max_value - 0.2 * magnitude;
         i += magnitude) {
      double tick = bar_ratio * (i - min_value) + ymin_bar;
      cairo_move_to(cr, xmax_bar - bar_width / 4, ly - tick);
      cairo_line_to(cr, xmax_bar, ly - tick);
      cairo_stroke(cr);
      cairo_move_to(
        cr,
        xmax_bar + bar_width / 4,
        ly - tick + (font_size / 2.0));
      cairo_show_text(cr, std::to_string(i).c_str());
    }

    // Negative ticks
    for (long long i = -magnitude; i > min_value + 0.2 * magnitude;
         i -= magnitude) {
      double tick = bar_ratio * (i - min_value) + ymin_bar;
      cairo_move_to(cr, xmax_bar - bar_width / 4, ly - tick);
      cairo_line_to(cr, xmax_bar, ly - tick);
      cairo_stroke(cr);
      cairo_move_to(
        cr,
        xmax_bar + bar_width / 4,
        ly - tick + (font_size / 2.0));
      cairo_show_text(cr, std::to_string(i).c_str());
    }
  }
}

// This function creates a simple SVG/PS file with a density bar
// void write_density_bar_image(std::string filename, const bool
// image_format_ps)
// {
//   // Create a cairo surface
//   cairo_surface_t *surface;

//   // Create a cairo surface
//   image_format_ps
//     ? surface = cairo_ps_surface_create(filename.c_str(), 80, 200)
//     : surface = cairo_svg_surface_create(filename.c_str(), 80, 200);
//   cairo_t *cr = cairo_create(surface);

//   // Write header
//   if (image_format_ps) {
//     write_ps_header(filename, surface);
//   }

//   write_density_bar_to_cairo_surface(
//     -3,
//     0,
//     3,
//     cr,
//     Bbox(20.0, 15.0, 35.0, 165.0),
//     200);

//   cairo_show_page(cr);
//   cairo_surface_destroy(surface);
//   cairo_destroy(cr);
// }

void InsetState::write_density_image(
  const std::string filename,
  const double *density,
  const bool plot_graticule_heatmap,
  const bool image_format_ps)
{
  // Whether to draw bar on the cairo surface
  const bool draw_bar = true;
  cairo_surface_t *surface;

  // Create a cairo surface
  image_format_ps
    ? surface = cairo_ps_surface_create(filename.c_str(), lx_, ly_)
    : surface = cairo_svg_surface_create(filename.c_str(), lx_, ly_);
  cairo_t *cr = cairo_create(surface);

  // Write header
  if (image_format_ps) {
    write_ps_header(filename, surface);
  }

  const Bbox bbox_bar = get_bbox_bar(15, 150);

  // double each_graticule_cell_area_km = graticule_cell_area_km(0, 0);
  cairo_set_line_width(cr, 0);

  // Determine range of densities
  double dens_min = dens_min_;
  double dens_mean = dens_mean_;
  double dens_max = dens_max_;

  // Crop it too
  if (plot_graticule_heatmap) {
    unsigned int cell_width = 1;

    // Clip to shape and print in sequential scale
    for (const auto &gd : geo_divs_) {
      for (const auto &pwh : gd.polygons_with_holes()) {
        const auto ext_ring = pwh.outer_boundary();
        cairo_new_path(cr);
        cairo_move_to(cr, ext_ring[0].x(), ly_ - ext_ring[0].y());

        // Plot each point in exterior ring
        for (unsigned int i = 1; i < ext_ring.size(); ++i) {
          cairo_line_to(cr, ext_ring[i].x(), ly_ - ext_ring[i].y());
        }

        // Close entire path
        cairo_close_path(cr);
        cairo_clip(cr);

        Bbox gd_bbox = gd.bbox();

        // Iterate over graticule cells
        for (unsigned int i = gd_bbox.xmin(); i < gd_bbox.xmax();
             i += cell_width) {
          for (unsigned int j = gd_bbox.ymin(); j < gd_bbox.ymax();
               j += cell_width) {

            // Values here used are "Max target area per km" and
            // "Min target area per km", which is obtained by running the
            // code with the "plot_graticule_heatmap" -h flag set to true
            Color color =
              graticule_cell_color(density[i * ly_ + j], dens_max, dens_min);

            // Get four points of the square
            double x_min = i - 0.5 * sq_overlap;
            double y_min = j - 0.5 * sq_overlap;
            double x_max = i + 1 + 0.5 * sq_overlap;
            double y_max = j + 1 + 0.5 * sq_overlap;

            cairo_move_to(cr, x_min, ly_ - y_min);
            cairo_line_to(cr, x_max, ly_ - y_min);
            cairo_line_to(cr, x_max, ly_ - y_max);
            cairo_line_to(cr, x_min, ly_ - y_max);

            cairo_set_source_rgb(cr, color.r, color.g, color.b);
            cairo_fill(cr);
            cairo_set_source_rgb(cr, 0, 0, 0);
            cairo_stroke(cr);

            // Fill the graticule polygon with color
            cairo_fill_preserve(cr);
            cairo_stroke(cr);
          }
        }

        // Remove GeoDiv clip
        cairo_reset_clip(cr);
      }
    }
  } else {
    for (unsigned int i = 0; i < lx_; ++i) {
      for (unsigned int j = 0; j < ly_; ++j) {
        Color color =
          heatmap_color(density[i * ly_ + j], dens_min, dens_mean, dens_max);

        // Get four points of the square
        double x_min = i - 0.5 * sq_overlap;
        double y_min = j - 0.5 * sq_overlap;
        double x_max = i + 1 + 0.5 * sq_overlap;
        double y_max = j + 1 + 0.5 * sq_overlap;

        cairo_move_to(cr, x_min, ly_ - y_min);
        cairo_line_to(cr, x_max, ly_ - y_min);
        cairo_line_to(cr, x_max, ly_ - y_max);
        cairo_line_to(cr, x_min, ly_ - y_max);

        cairo_set_source_rgb(cr, color.r, color.g, color.b);
        cairo_fill(cr);
        cairo_set_source_rgb(cr, 0, 0, 0);
        cairo_stroke(cr);
      }
    }
  }
  write_polygons_to_cairo_surface(cr, false, false, false);

  if (draw_bar && !plot_graticule_heatmap) {

    write_density_bar_to_cairo_surface(
      dens_min - dens_mean,
      0,
      dens_max - dens_mean,
      cr,
      bbox_bar,
      ly_);

    std::string bar_filename =
      "bar_" + filename.substr(0, filename.size() - 2) + "svg";
    // std::string bar_filename = "bar_" + filename;

    // Create a cairo bar_surface
    cairo_surface_t *bar_surface =
      cairo_svg_surface_create(bar_filename.c_str(), 160, 400);
    // cairo_surface_t *bar_surface;
    // image_format_ps
    //   ? bar_surface = cairo_ps_surface_create(bar_filename.c_str(), 160, 400)
    //   : bar_surface = cairo_svg_surface_create(bar_filename.c_str(), 160,
    //   400);
    cairo_t *bar_cr = cairo_create(bar_surface);

    // Write header
    // if (image_format_ps) {
    //   write_ps_header(bar_filename, bar_surface);
    // }

    write_density_bar_to_cairo_surface(
      dens_min - dens_mean,
      0,
      dens_max - dens_mean,
      bar_cr,
      Bbox(75.0, 200.0, 95.0, 350.0),
      ly_);

    cairo_show_page(bar_cr);
    cairo_surface_destroy(bar_surface);
    cairo_destroy(bar_cr);
  }

  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

void InsetState::write_intersections_image(
  unsigned int res,
  const bool image_format_ps)
{
  std::string filename = inset_name() + "_intersections_" +
                         std::to_string(n_finished_integrations());

  // Update extension
  image_format_ps ? filename += ".ps" : filename += ".svg";

  // Calculating intersections
  std::vector<Segment> intersections = intersecting_segments(res);
  cairo_surface_t *surface;

  // Create a cairo surface
  image_format_ps
    ? surface = cairo_ps_surface_create(filename.c_str(), lx_, ly_)
    : surface = cairo_svg_surface_create(filename.c_str(), lx_, ly_);
  cairo_t *cr = cairo_create(surface);

  // Write header
  if (image_format_ps) {
    write_ps_header(filename, surface);
  }

  write_polygons_to_cairo_surface(cr, false, false, false);

  cairo_set_line_width(cr, 0.0001 * std::min(lx_, ly_));

  for (auto seg : intersections) {
    // Move to starting coordinates
    cairo_move_to(cr, seg[0][0], ly_ - seg[0][1]);

    // Draw line
    cairo_line_to(cr, seg[1][0], ly_ - seg[1][1]);

    // line with red and stroke
    cairo_set_source_rgb(cr, 1, 0, 0);
    cairo_stroke(cr);
  }
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}
