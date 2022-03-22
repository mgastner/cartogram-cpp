#include "constants.h"
#include "colors.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include <cairo/cairo.h>
#include <cairo/cairo-ps.h>
#include <cairo/cairo-svg.h>

void write_ps_header(const std::string filename,
                     cairo_surface_t *surface)
{
  const std::string title = "%%Title: " + filename;
  cairo_ps_surface_dsc_comment(surface, title.c_str());
  cairo_ps_surface_dsc_comment(surface,
                               "%%Creator: Michael T. Gastner et al.");
  cairo_ps_surface_dsc_comment(surface,
                               "%%For: Humanity");
  cairo_ps_surface_dsc_comment(surface,
                               "%%Copyright: License CC BY");
  cairo_ps_surface_dsc_comment(surface,
                               "%%Magnification: 1.0000");
}

// TODO: IS THERE A CGAL WAY OF DETERMINING WHETHER THE LABEL'S BOUNDING
//       BOX IS COMPLETELY CONTAINED IN THE POLYGON?

// TODO: SHOULD THE CRITERION FOR PRINTING A LABEL BE THAT IT FITS INSIDE THE
//       POLYGON WITH HOLES? THAT CRITERION WOULD BE MORE RESTRICTIVE THAN
//       FITTING INSIDE THE EXTERIOR RING.
bool all_points_inside_exterior_ring(const std::vector<Point> pts,
                                     const Polygon_with_holes pwh)
{
  for (const auto &pt : pts)
  {
    if (pwh.outer_boundary().has_on_unbounded_side(pt))
    {
      return false;
    }
  }
  return true;
}

double font_size(cairo_t *cr,
                 const char *label,
                 const Point label_pt,
                 const GeoDiv gd)
{
  for (double fsize = max_font_size; fsize >= min_font_size; fsize -= 0.5)
  {
    cairo_set_font_size(cr, fsize);
    cairo_text_extents_t extents;
    cairo_text_extents(cr, label, &extents);
    const Polygon_with_holes largest_pwh = gd.largest_polygon_with_holes();

    // Bounding box of the label
    const CGAL::Bbox_2 bb(label_pt.x() - 0.5 * extents.width,
                          label_pt.y() - 0.5 * extents.height,
                          label_pt.x() + 0.5 * extents.width,
                          label_pt.y() + 0.5 * extents.height);

    // Vector of bounding-box edge points
    std::vector<Point> bb_edge_points;
    for (unsigned int i = 0; i <= 1; ++i)
    {
      for (unsigned int j = 0; j <= 5; ++j)
      {
        bb_edge_points.push_back(
            Point((j * bb.xmin() + (5 - j) * bb.xmax()) / 5,
                  (i * bb.ymin() + (1 - i) * bb.ymax())));
      }
    }
    if (all_points_inside_exterior_ring(bb_edge_points, largest_pwh))
    {
      return fsize;
    }
  }
  return 0.0;
}

void write_labels_to_cairo_surface(cairo_t *cr,
                                   const InsetState *inset_state)
{
  const unsigned int ly = inset_state->ly();
  for (const auto &gd : inset_state->geo_divs())
  {
    const std::string label = inset_state->label_at(gd.id());
    const char *const label_char = label.c_str();
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_select_font_face(cr,
                           "sans-serif",
                           CAIRO_FONT_SLANT_NORMAL,
                           CAIRO_FONT_WEIGHT_NORMAL);
    const Point label_pt = gd.point_on_surface_of_geodiv();

    // Get size of label
    const double fsize = font_size(cr, label_char, label_pt, gd);
    cairo_text_extents_t extents;

    // Draw label only if appropriate size is found
    if (fsize > 0.0)
    {
      cairo_set_font_size(cr, fsize);
      cairo_text_extents(cr, label_char, &extents);
      const double x =
          label_pt.x() - (extents.width / 2 + extents.x_bearing);
      const double y =
          ly - label_pt.y() - (extents.height / 2 + extents.y_bearing);
      cairo_move_to(cr, x, y);
      cairo_show_text(cr, label_char);
    }
  }
}

// Returns edge points of a graticule cell and treats
// them as a polygon.
Polygon graticule_cell_edge_points(unsigned int x,
                                   unsigned int y,
                                   InsetState *inset_state,
                                   unsigned int cell_width = graticule_width,
                                   bool plot_equal_area_map = false)
{
  Polygon cell_edge_points;
  boost::multi_array<XYPoint, 2> &cum_proj = plot_equal_area_map ? *inset_state->ref_to_original_proj() : *inset_state->ref_to_cum_proj();

  // Horizontal lower edge points
  for (unsigned int i = x; i < x + cell_width; ++i)
  {
    double x_coor_trans = cum_proj[i][y].x;
    double y_coor_trans = cum_proj[i][y].y;
    cell_edge_points.push_back(Point(x_coor_trans, y_coor_trans));
  }

  // Vertical right edge points
  for (unsigned int i = y; i < y + cell_width; ++i)
  {
    double x_coor_trans = cum_proj[x + cell_width][i].x;
    double y_coor_trans = cum_proj[x + cell_width][i].y;
    cell_edge_points.push_back(Point(x_coor_trans, y_coor_trans));
  }

  // Horizontal upper edge points
  for (unsigned int i = x + cell_width; i > x; --i)
  {
    double x_coor_trans = cum_proj[i][y + cell_width].x;
    double y_coor_trans = cum_proj[i][y + cell_width].y;
    cell_edge_points.push_back(Point(x_coor_trans, y_coor_trans));
  }

  // Vertical left edge points
  for (unsigned int i = y + cell_width; i > y; --i)
  {
    double x_coor_trans = cum_proj[x][i].x;
    double y_coor_trans = cum_proj[x][i].y;
    cell_edge_points.push_back(Point(x_coor_trans, y_coor_trans));
  }

  // Complete the polygon by making first and last point same
  cell_edge_points.push_back(cell_edge_points[0]);

  return cell_edge_points;
}

// Returns graticule cell area based on edge points
double graticule_cell_area(unsigned int x,
                           unsigned int y,
                           InsetState *inset_state,
                           unsigned int cell_width)
{

  // Taking absolule to ensure we get the area irrespective of direction
  return abs(graticule_cell_edge_points(x, y, inset_state, cell_width).area());
}

// Returns the largest and smallest graticule cell area to be used for
// graticule heatmap generation
std::pair<double, double> max_and_min_graticule_cell_area(InsetState *inset_state,
                                                          unsigned int cell_width)
{
  unsigned int lx = inset_state->lx();
  unsigned int ly = inset_state->ly();

  // Initialize max and min area
  double max_area = -dbl_inf;
  double min_area = dbl_inf;

  // Iterate over graticule cells
  for (unsigned int i = 0; i < lx - cell_width; i += cell_width)
  {
    for (unsigned int j = 0; j < ly - cell_width; j += cell_width)
    {
      const double area = graticule_cell_area(i, j, inset_state, cell_width);
      max_area = std::max(max_area, area);
      min_area = std::min(min_area, area);
    }
  }
  return std::make_pair(max_area, min_area);
}

std::pair<Point, Point> max_and_min_graticule_cell_area_index(InsetState *inset_state,
                                                              unsigned int cell_width)
{
  unsigned int lx = inset_state->lx();
  unsigned int ly = inset_state->ly();

  // Initialize max and min area
  double max_area = -dbl_inf;
  double min_area = dbl_inf;
  unsigned int max_i = 0;
  unsigned int max_j = 0;
  unsigned int min_i = 0;
  unsigned int min_j = 0;

  // Iterate over graticule cells
  for (unsigned int i = 0; i < lx - cell_width; i += cell_width)
  {
    for (unsigned int j = 0; j < ly - cell_width; j += cell_width)
    {
      const double area = graticule_cell_area(i, j, inset_state, cell_width);
      if (area > max_area)
      {
        max_area = area;
        max_i = i;
        max_j = j;
      }
      if (area < min_area)
      {
        min_area = area;
        min_i = i;
        min_j = j;
      }
    }
  }
  return std::make_pair(Point(max_i, max_j), Point(min_i, min_j));
}

void graticule_cell_color(const double area,
                          const double max_area,
                          const double min_area,
                          double *r,
                          double *g,
                          double *b)
{
  // Assign possible categories for red, green, blue
  const double red[] = {
      1.000, 0.996, 0.988, 0.988, 0.984, 0.937, 0.796, 0.647, 0.404};
  const double green[] = {
      0.961, 0.878, 0.733, 0.572, 0.416, 0.231, 0.094, 0.058, 0.000};
  const double blue[] = {
      0.941, 0.824, 0.631, 0.447, 0.290, 0.173, 0.114, 0.082, 0.050};

  // Normalize area to [0,1] and make it logarithmic
  double ratio = (log(area) - log(min_area)) / (log(max_area) - log(min_area));

  // Determine color category
  double category = fmax(0, ratio * 8 - 10e-6);
  double xmin = floor(category);
  double xmax = ceil(category);

  if (area == max_area)
  {
    *r = red[8];
    *g = green[8];
    *b = blue[8];
    return;
  }
  else if (area == min_area)
  {
    *r = red[0];
    *g = green[0];
    *b = blue[0];
    return;
  }
  else
  {
    *r = red[int(xmin)] + (red[int(xmax)] - red[int(xmin)]) * (category - xmin);
    *g = green[int(xmin)] + (green[int(xmax)] - green[int(xmin)]) * (category - xmin);
    *b = blue[int(xmin)] + (blue[int(xmax)] - blue[int(xmin)]) * (category - xmin);
  }
}

std::vector<std::vector<Color_dbl>> graticule_cell_colors(InsetState *inset_state,
                                                          unsigned int cell_width)
{
  unsigned int lx = inset_state->lx();
  unsigned int ly = inset_state->ly();

  // Initialize max and min area
  double max_area, min_area;
  std::tie(max_area, min_area) = max_and_min_graticule_cell_area(inset_state,
                                                                 cell_width);

  // Initialize colors
  std::vector<std::vector<Color_dbl>> colors(lx - cell_width,
                                             std::vector<Color_dbl>(ly - cell_width));

  // Iterate over graticule cells
  for (unsigned int i = 0; i < lx - cell_width; i += cell_width)
  {
    for (unsigned int j = 0; j < ly - cell_width; j += cell_width)
    {
      const double area = graticule_cell_area(i, j, inset_state, cell_width);
      double r, g, b;
      graticule_cell_color(area, max_area, min_area, &r, &g, &b);
      Color_dbl color = {r, g, b};
      colors[i][j] = color;
    }
  }
  return colors;
}

void write_graticule_heatmap_bar_to_cairo_surface(double min_value,
                                                  double max_value,
                                                  cairo_t *cr,
                                                  Bbox bbox_bar,
                                                  std::vector<std::pair<double, double>> major_ticks,
                                                  std::vector<std::pair<double, double>> minor_ticks,
                                                  const unsigned int ly)
{
  const int n_gradident_bars = 500;

  // get bar coordinates
  const double xmin_bar = bbox_bar.xmin();
  const double xmax_bar = bbox_bar.xmax();
  const double ymin_bar = bbox_bar.ymin();
  const double ymax_bar = bbox_bar.ymax();

  const double bar_width = xmax_bar - xmin_bar;

  // calculate individual bar gradient segment property
  const double gradient_segment_height = (ymax_bar - ymin_bar) / n_gradident_bars;
  const double gradient_segment_value = (max_value - min_value) / n_gradident_bars;

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

  for (double y = ymin_bar; y <= ymax_bar; y += gradient_segment_height)
  {
    double r, g, b;
    graticule_cell_color(exp(value_at_gradient_segment), exp(max_value),
                         exp(min_value), &r, &g, &b);
    cairo_set_source_rgb(cr, r, g, b);
    cairo_rectangle(cr, xmin_bar, ly - y, bar_width, gradient_segment_height);
    cairo_fill(cr);
    value_at_gradient_segment += gradient_segment_value;
  }

  // Set font properties
  cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size(cr, 8);

  // Draw the ticks and nice_numbers
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_set_line_width(cr, .9);
  for (auto tick : major_ticks)
  {
    double area = tick.first;
    int NiceNumber = tick.second;
    if (area > min_value and area < max_value)
    {
      double y = ((log(area) - log(min_value)) / (log(max_value) - log(min_value))) * (ymax_bar - ymin_bar) + ymin_bar;
      cairo_move_to(cr, xmax_bar - 6, ly - y);
      cairo_line_to(cr, xmax_bar, ly - y);
      cairo_move_to(cr, xmax_bar + 1, ly - y);
      cairo_show_text(cr, std::to_string(NiceNumber).c_str());
      cairo_stroke(cr);
    }
  }
  cairo_set_line_width(cr, .7);
  for (auto ticks : minor_ticks)
  {
    double area = ticks.first;
    if (area > min_value and area < max_value)
    {
      double y = ((log(area) - log(min_value)) / (log(max_value) - log(min_value))) *
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
  cairo_move_to(cr, (xmin_bar + xmax_bar) / 2 - 16, ly - ymax_bar - 13.5);
  cairo_show_text(cr, bar_text_top.c_str());
  cairo_move_to(cr, (xmin_bar + xmax_bar) / 2 - 5.5, ly - ymax_bar - 5);
  cairo_show_text(cr, bar_text_bottom.c_str());
}

void trim_graticule_heatmap(cairo_t *cr,
                            double padding,
                            InsetState *inset_state)
{
  // Canvas dimension
  double lx = (double)inset_state->lx();
  double ly = (double)inset_state->ly();
  Bbox bbox = inset_state->bbox();
  double xmin = bbox.xmin() - padding;
  double xmax = bbox.xmax() + padding;
  double ymin = bbox.ymin() - padding;
  double ymax = bbox.ymax() + padding;

  // Color white outside bbox
  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_rectangle(cr, 0, 0, xmin, ymax);
  cairo_rectangle(cr, xmax, 0, lx - xmax, ymax);
  cairo_rectangle(cr, 0, ymax, lx, ly - ymax);
  cairo_rectangle(cr, 0, 0, lx, ymin);
  cairo_fill(cr);
}

// Writes graticule cells and colors them if required
void write_graticules_to_cairo_surface(cairo_t *cr,
                                       InsetState *inset_state)
{
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();

  // Set line width of graticule lines
  cairo_set_line_width(cr, 5e-4 * std::min(lx, ly));
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);

  // Iterate over graticule cells
  for (unsigned int i = 0; i < lx - graticule_width; i += graticule_width)
  {
    for (unsigned int j = 0; j < ly - graticule_width; j += graticule_width)
    {

      // Draw graticule cell by connecting edge points
      const Polygon cell_edge_points = graticule_cell_edge_points(i, j, inset_state);
      cairo_move_to(cr, cell_edge_points[0].x(), ly - cell_edge_points[0].y());
      for (unsigned int k = 1; k < cell_edge_points.size(); ++k)
      {
        cairo_line_to(cr, cell_edge_points[k].x(), ly - cell_edge_points[k].y());
      }
      cairo_stroke(cr);
    }
  }
}

void write_graticule_colors_to_cairo_surface(cairo_t *cr,
                                             InsetState *inset_state,
                                             bool plot_equal_area_map)
{
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  unsigned int cell_width = 1;

  // Get colors
  const auto colors = graticule_cell_colors(inset_state, cell_width);

  // Set line width of graticule lines
  cairo_set_line_width(cr, 5e-6 * std::min(lx, ly));

  // Iterate over graticule cells
  for (unsigned int i = 0; i < lx - cell_width; i += cell_width)
  {
    for (unsigned int j = 0; j < ly - cell_width; j += cell_width)
    {

      // Set color of the border of the graticule polygon
      cairo_set_source_rgb(cr, colors[i][j].r, colors[i][j].g, colors[i][j].b);

      // Draw graticule cell by connecting edge points
      const Polygon cell_edge_points = graticule_cell_edge_points(i, j, inset_state, cell_width,
                                                               plot_equal_area_map);
      cairo_move_to(cr, cell_edge_points[0].x(), ly - cell_edge_points[0].y());
      for (unsigned int k = 1; k < cell_edge_points.size(); ++k)
      {
        cairo_line_to(cr, cell_edge_points[k].x(), ly - cell_edge_points[k].y());
      }

      // Fill the graticule polygon with color
      cairo_fill_preserve(cr);
      cairo_stroke(cr);
    }
  }
}

void write_polygons_to_cairo_surface(cairo_t *cr,
                                     const bool fill_polygons,
                                     const bool colors,
                                     const bool plot_equal_area_map,
                                     InsetState *inset_state)
{
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  cairo_set_line_width(cr, 1e-3 * std::min(lx, ly));

  // Draw cartogram polygons or equal area map polygons
  std::vector<GeoDiv> geo_divs = plot_equal_area_map ? inset_state->geo_divs_original() : inset_state->geo_divs();

  // Draw the shapes
  for (const auto &gd : geo_divs)
  {
    for (const auto &pwh : gd.polygons_with_holes())
    {
      const Polygon ext_ring = pwh.outer_boundary();
      cairo_move_to(cr, ext_ring[0].x(), ly - ext_ring[0].y());

      // Plot each point in exterior ring
      for (unsigned int i = 1; i < ext_ring.size(); ++i)
      {
        cairo_line_to(cr, ext_ring[i].x(), ly - ext_ring[i].y());
      }

      // Close the exterior ring
      cairo_close_path(cr);

      // Draw holes
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h)
      {
        cairo_move_to(cr, (*h)[0].x(), ly - (*h)[0].y());
        const size_t hsize = (*h).size();
        for (unsigned int i = 1; i <= hsize; ++i)
        {
          cairo_line_to(cr, (*h)[i % hsize].x(), ly - (*h)[i % hsize].y());
        }
      }
      if (colors || fill_polygons)
      {
        if (inset_state->is_input_target_area_missing(gd.id()))
        {

          // Fill path with dark gray
          cairo_set_source_rgb(cr, 0.9375, 0.9375, 0.9375);
        }
        else if (colors)
        {

          // Get color
          const Color col = inset_state->color_at(gd.id());

          // Fill path
          cairo_set_source_rgb(cr,
                               col.r / 255.0,
                               col.g / 255.0,
                               col.b / 255.0);
        }
        else if (fill_polygons)
        {

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
void write_map_image(const std::string filename,
                     const bool fill_polygons,
                     const bool plot_graticule,
                     const bool image_format_ps,
                     InsetState *inset_state)
{
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();

  // Check whether the map has all GeoDivs colored
  const bool colors =
      (inset_state->colors_size() == inset_state->n_geo_divs());
  cairo_surface_t *surface;

  // Create a cairo surface
  image_format_ps ? surface = cairo_ps_surface_create(filename.c_str(), lx, ly) : surface = cairo_svg_surface_create(filename.c_str(), lx, ly);
  cairo_t *cr = cairo_create(surface);

  // Write header
  if (image_format_ps)
  {

    write_ps_header(filename, surface);
  }

  // Draw polygons with color
  write_polygons_to_cairo_surface(cr,
                                  fill_polygons,
                                  colors,
                                  false,
                                  inset_state);

  // Place labels
  write_labels_to_cairo_surface(cr, inset_state);

  // Draw graticule without color
  if (plot_graticule)
  {
    write_graticules_to_cairo_surface(cr, inset_state);
  }
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

Polygon transform_to_albers_coor(Polygon edge_points,
                                 const InsetState *inset_state)
{
  const double latt_const = inset_state->latt_const();

  Transformation scale(CGAL::SCALING, latt_const);

  Polygon cell_edge_points_albers = transform(scale, edge_points);

  return cell_edge_points_albers;
}

double albers_area_to_earth_area(const double albers_area)
{

  return (albers_area * earth_surface_area) / (4 * pi);
}

double graticule_cell_area_km(const unsigned int i,
                              const unsigned int j,
                              InsetState *inset_state)
{

  Polygon cell_edge_points = graticule_cell_edge_points(i, j, inset_state,
                                                        1,
                                                        true);
  const Polygon cell_edge_points_albers = transform_to_albers_coor(cell_edge_points,
                                                                   inset_state);
  const double cell_area = cell_edge_points_albers.area();

  const double cell_area_km = albers_area_to_earth_area(cell_area);

  return cell_area_km;
}

double graticule_cell_target_area(const unsigned int i,
                                  const unsigned int j,
                                  const double total_target_area,
                                  const double total_inset_area,
                                  InsetState *inset_state)
{
  const Polygon cell_edge_points = graticule_cell_edge_points(i, j, inset_state,
                                                              1,
                                                              false);
  const double cell_area = cell_edge_points.area();

  const double cell_target_area = (cell_area * total_target_area) / total_inset_area;

  return cell_target_area;
}

double graticule_cell_target_area_per_km(const unsigned int i,
                                         const unsigned int j,
                                         const double total_target_area,
                                         const double total_inset_area,
                                         InsetState *inset_state)
{
  const double cell_target_area = graticule_cell_target_area(i, j,
                                                             total_target_area,
                                                             total_inset_area,
                                                             inset_state);
  const double cell_area_km = graticule_cell_area_km(i, j, inset_state);

  const double cell_target_area_per_km = cell_target_area / cell_area_km;

  return cell_target_area_per_km;
}

std::vector<std::pair<double, double>> get_major_ticks(const double min_target_area_per_km,
                                                       const double max_target_area_per_km,
                                                       const double min_area_cell_point_area,
                                                       const double max_area_cell_point_area,
                                                       std::vector<int> nice_numbers)

{
  std::vector<std::pair<double, double>> ticks;
  for (auto niceNumber : nice_numbers)
  {
    double NiceNumberRatio = (niceNumber - min_target_area_per_km) / (max_target_area_per_km - min_target_area_per_km);
    double area = min_area_cell_point_area + NiceNumberRatio *
                                                 (max_area_cell_point_area - min_area_cell_point_area);
    ticks.push_back(std::make_pair(area, niceNumber));
  }
  return ticks;
}

std::vector<std::pair<double, double>> get_minor_ticks(int n_ticks_per_major,
                                                       const double min_target_area_per_km,
                                                       const double max_target_area_per_km,
                                                       const double min_area_cell_point_area,
                                                       const double max_area_cell_point_area,
                                                       std::vector<int> nice_numbers)

{
  n_ticks_per_major += 2;
  std::vector<std::pair<double, double>> minor_ticks;
  const int n_major_ticks = nice_numbers.size();
  for (int i = 0; i < n_major_ticks - 1; i++)
  {
    double first_major_tick = (double)nice_numbers[i];
    double second_major_tick = (double)nice_numbers[i + 1];
    double minor_tick_ratio = (second_major_tick - first_major_tick) / (n_ticks_per_major - 1);
    for (int j = 1; j < n_ticks_per_major - 1; j++)
    {
      double minor_tick = first_major_tick + minor_tick_ratio * j;
      double NiceNumberRatio = (minor_tick - min_target_area_per_km) /
                               (max_target_area_per_km - min_target_area_per_km);
      double area = min_area_cell_point_area + NiceNumberRatio *
                                                   (max_area_cell_point_area - min_area_cell_point_area);
      minor_ticks.push_back(std::make_pair(area, minor_tick));
    }
  }
  return minor_ticks;
}

std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>>
get_ticks(const int n_ticks_per_major,
          const double min_target_area_per_km,
          const double max_target_area_per_km,
          const double min_area_cell_point_area,
          const double max_area_cell_point_area,
          std::vector<int> nice_numbers)
{
  std::vector<std::pair<double, double>> major_ticks = get_major_ticks(min_target_area_per_km,
                                                                       max_target_area_per_km,
                                                                       min_area_cell_point_area,
                                                                       max_area_cell_point_area,
                                                                       nice_numbers);
  std::vector<std::pair<double, double>> minor_ticks = get_minor_ticks(n_ticks_per_major,
                                                                       min_target_area_per_km,
                                                                       max_target_area_per_km,
                                                                       min_area_cell_point_area,
                                                                       max_area_cell_point_area,
                                                                       nice_numbers);
  std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>>
      ticks(major_ticks, minor_ticks);
  return ticks;
}

Bbox get_bbox_bar(const double bar_width,
                  const double bar_height,
                  InsetState *inset_state)
{

  const Bbox bbox = inset_state->bbox();

  // Position the bar 25 pixels to the right of the bbox
  const double xmin_bar = bbox.xmax() + 35;
  const double xmax_bar = xmin_bar + bar_width;

  // Position the bar at the middle of the bbox y coordinates
  double ymid_bar = (bbox.ymax() + bbox.ymin()) / 2 - 25;
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
  while (NiceNumber < max_target_area_per_km)
  {
    NiceNumber = NiceNumber * 10;
    nice_numbers.push_back(NiceNumber);
  }
  return nice_numbers;
}

// Outputs a SVG/PS file of graticule heatmap
void write_graticule_heatmap_image(const std::string filename,
                                   const bool plot_equal_area_map,
                                   const bool image_format_ps,
                                   InsetState *inset_state)
{
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  cairo_surface_t *surface;

  // Create a cairo surface
  image_format_ps ? surface = cairo_ps_surface_create(filename.c_str(), lx, ly) : surface = cairo_svg_surface_create(filename.c_str(), lx, ly);
  cairo_t *cr = cairo_create(surface);

  // Get inset areas
  const double total_target_area = inset_state->total_target_area();
  const double total_inset_area = inset_state->total_inset_area();

  const Bbox bbox_bar = get_bbox_bar(15, 150, inset_state);

  // Get the max and min graticule cell area points
  Point max_area_cell_point, min_area_cell_point;

  std::tie(max_area_cell_point, min_area_cell_point) = max_and_min_graticule_cell_area_index(inset_state, 1);

  const double max_area_cell_point_area = graticule_cell_area(max_area_cell_point.x(),
                                                              max_area_cell_point.y(),
                                                              inset_state, 1);
  const double min_area_cell_point_area = graticule_cell_area(min_area_cell_point.x(),
                                                              min_area_cell_point.y(),
                                                              inset_state, 1);

  // Get the max and min graticule cell target area per km
  double max_target_area_per_km = graticule_cell_target_area_per_km(max_area_cell_point.x(),
                                                                    max_area_cell_point.y(),
                                                                    total_target_area,
                                                                    total_inset_area,
                                                                    inset_state);
  double min_target_area_per_km = graticule_cell_target_area_per_km(min_area_cell_point.x(),
                                                                    min_area_cell_point.y(),
                                                                    total_target_area,
                                                                    total_inset_area,
                                                                    inset_state);
  std::cerr << std::endl;
  std::cerr << "Max target area per km: " << max_target_area_per_km << std::endl;
  std::cerr << "Min target area per km: " << min_target_area_per_km << std::endl;

  std::vector<int> nice_numbers = get_nice_numbers_for_bar(max_target_area_per_km);

  std::vector<std::pair<double, double>> major_ticks, minor_ticks;

  std::tie(major_ticks, minor_ticks) = get_ticks(10, min_target_area_per_km, max_target_area_per_km,
                                                 min_area_cell_point_area, max_area_cell_point_area,
                                                 nice_numbers);

  // Write header
  if (image_format_ps)
  {
    write_ps_header(filename, surface);
  }

  // Draw colors
  write_graticule_colors_to_cairo_surface(cr, inset_state, plot_equal_area_map);

  // Draw polygons without color
  write_polygons_to_cairo_surface(cr,
                                  false,
                                  false,
                                  plot_equal_area_map,
                                  inset_state);

  trim_graticule_heatmap(cr, 20, inset_state);

  write_graticule_heatmap_bar_to_cairo_surface(min_area_cell_point_area,
                                               max_area_cell_point_area,
                                               cr, bbox_bar, major_ticks,
                                               minor_ticks, ly);
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

// Functions to show a scalar field called "density" as a heat map
double interpolate_for_heatmap(const double x,
                               const double xmin,
                               const double xmax,
                               const double ymin,
                               const double ymax)
{
  return ((x - xmin) * ymax + (xmax - x) * ymin) / (xmax - xmin);
}

void heatmap_color(const double dens,
                   const double dens_min,
                   const double dens_mean,
                   const double dens_max,
                   double *r,
                   double *g,
                   double *b)
{
  // Assign possible categories for red, green, blue
  const double red[] = {
      0.33, 0.55, 0.75, 0.87, 0.96, 0.99, 0.78, 0.50, 0.21, 0.00, 0.00};
  const double green[] = {
      0.19, 0.32, 0.51, 0.76, 0.91, 0.96, 0.92, 0.80, 0.59, 0.40, 0.24};
  const double blue[] = {
      0.02, 0.04, 0.18, 0.49, 0.76, 0.89, 0.90, 0.76, 0.56, 0.37, 0.19};
  double xmin, xmax;
  int color_category;

  // Choose color category
  if (dens > dens_max)
  {
    *r = red[0];
    *g = green[0];
    *b = blue[0];
    return;
  }
  else if (dens > dens_mean)
  {
    color_category = 5 * (dens_max - dens) / (dens_max - dens_mean);
    xmax = dens_max - 0.2 * color_category * (dens_max - dens_mean);
    xmin = xmax - 0.2 * (dens_max - dens_mean);

    // Assign color category 0 if dens_max and dens are very close
    color_category = std::max(color_category, 0);
  }
  else if (dens > dens_min)
  {
    color_category = 5 * (dens_mean - dens) / (dens_mean - dens_min) + 5;
    xmax = dens_mean - 0.2 * (color_category - 5) * (dens_mean - dens_min);
    xmin = xmax - 0.2 * (dens_mean - dens_min);

    // Assign color category 9 if dens_min and dens are very close
    color_category = std::min(color_category, 9);
  }
  else
  {
    *r = red[10];
    *g = green[10];
    *b = blue[10];
    return;
  }
  *r = interpolate_for_heatmap(dens,
                               xmin,
                               xmax,
                               red[color_category + 1],
                               red[color_category]);
  *g = interpolate_for_heatmap(dens,
                               xmin,
                               xmax,
                               green[color_category + 1],
                               green[color_category]);
  *b = interpolate_for_heatmap(dens,
                               xmin,
                               xmax,
                               blue[color_category + 1],
                               blue[color_category]);
}

// Function to show the density bar on the cairo surface
void write_density_bar_to_cairo_surface(const double min_value,
                                        const double mean_value,
                                        const double max_value,
                                        cairo_t *cr,
                                        Bbox bbox_bar,
                                        const unsigned int ly)
{
  const int n_gradident_bars = 500;

  // get bar coordinates
  const double xmin_bar = bbox_bar.xmin();
  const double xmax_bar = bbox_bar.xmax();
  const double ymin_bar = bbox_bar.ymin();
  const double ymax_bar = bbox_bar.ymax();

  const double bar_width = xmax_bar - xmin_bar;

  // position of mean line along bar
  const double ymean_bar = ((ymax_bar - ymin_bar) / (max_value - min_value)) * (mean_value - min_value) + ymin_bar;

  // calculate individual bar gradient segment property
  const double gradient_segment_height = (ymax_bar - ymin_bar) / n_gradident_bars;
  const double gradient_segment_value = (max_value - min_value) / n_gradident_bars;

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

  for (double y = ymin_bar; y <= ymax_bar; y += gradient_segment_height)
  {
    double r, g, b;
    heatmap_color(value_at_gradient_segment, min_value, mean_value, max_value, &r, &g, &b);
    cairo_set_source_rgb(cr, r, g, b);
    cairo_rectangle(cr, xmin_bar, ly - y, bar_width, gradient_segment_height);
    cairo_fill(cr);
    value_at_gradient_segment += gradient_segment_value;
  }

  // Draw the mean line
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_move_to(cr, xmin_bar + bar_width / 2, ly - ymean_bar);
  cairo_line_to(cr, xmax_bar, ly - ymean_bar);
  cairo_stroke(cr);

  // Set font properties
  cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size(cr, 10);

  // Write "Density" on top of bar
  cairo_move_to(cr, xmin_bar - bar_width / 2 - 1, ly - (ymax_bar + 15));
  cairo_show_text(cr, "Density");

  // Write "High" top right of bar
  cairo_move_to(cr, xmax_bar + 5, ly - ymax_bar + 3);
  cairo_show_text(cr, "High");

  // Write "Low" bottom right of bar
  cairo_move_to(cr, xmax_bar + 5, ly - ymin_bar + 3);
  cairo_show_text(cr, "Low");

  // Write "Mean" beside ymean_bar
  cairo_move_to(cr, xmax_bar + 5, ly - ymean_bar + 3);
  cairo_show_text(cr, "Mean");
}

// This function creates a simple SVG/PS file with a density bar
void write_density_bar_image(std::string filename,
                             const bool image_format_ps)
{

  // Create a cairo surface
  cairo_surface_t *surface;

  // Create a cairo surface
  image_format_ps ? surface = cairo_ps_surface_create(filename.c_str(), 80, 200) : surface = cairo_svg_surface_create(filename.c_str(), 80, 200);
  cairo_t *cr = cairo_create(surface);

  // Write header
  if (image_format_ps)
  {
    write_ps_header(filename, surface);
  }

  write_density_bar_to_cairo_surface(0,
                                     50,
                                     100,
                                     cr,
                                     Bbox(20.0, 15.0, 35.0, 165.0),
                                     200);

  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

void write_density_image(const std::string filename,
                         const double *density,
                         const bool plot_graticule_heatmap,
                         const bool image_format_ps,
                         InsetState *inset_state)
{
  // Whether to draw bar on the cairo surface
  const bool draw_bar = false;
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();

  cairo_surface_t *surface;

  // Create a cairo surface
  image_format_ps ? surface = cairo_ps_surface_create(filename.c_str(), lx, ly) : surface = cairo_svg_surface_create(filename.c_str(), lx, ly);
  cairo_t *cr = cairo_create(surface);

  // Write header
  if (image_format_ps)
  {

    write_ps_header(filename, surface);
  }

  const Bbox bbox_bar = get_bbox_bar(15, 150, inset_state);

  double each_graticule_cell_area_km = graticule_cell_area_km(0,
                                                              0,
                                                              inset_state);

  cairo_set_line_width(cr, 0);
  // Determine range of densities
  double dens_min = dbl_inf;
  double dens_mean = 0.0;
  double dens_max = -dbl_inf;
  const unsigned int n_grid_cells = inset_state->lx() * inset_state->ly();
  for (unsigned int k = 0; k < n_grid_cells; ++k)
  {
    dens_min = std::min(density[k], dens_min);
    dens_mean += density[k];
    dens_max = std::max(density[k], dens_max);
  }
  dens_mean /= n_grid_cells;

  for (unsigned int i = 0; i < inset_state->lx(); ++i)
  {
    for (unsigned int j = 0; j < inset_state->ly(); ++j)
    {
      double r, g, b;
      if (plot_graticule_heatmap)
      {
        double target_area_per_cell_km = density[i * inset_state->ly() + j] /
                                         each_graticule_cell_area_km;

        // Values here used are "Max target area per km" and
        // "Min target area per km", which is obtained by running the
        // code with the "plot_graticule_heatmap" -h flag set to true
        graticule_cell_color(target_area_per_cell_km,
                             642.872,
                             0.656938,
                             &r,
                             &g,
                             &b);
      }
      else
      {
        heatmap_color(density[i * inset_state->ly() + j],
                      dens_min,
                      dens_mean,
                      dens_max,
                      &r, &g, &b);
      }

      // Get four points of the square
      double x_min = i - 0.5 * sq_overlap;
      double y_min = j - 0.5 * sq_overlap;
      double x_max = x_min + 1.2;
      double y_max = y_min + 1.2;

      cairo_move_to(cr, x_min, ly - y_min);
      cairo_line_to(cr, x_max, ly - y_min);
      cairo_line_to(cr, x_max, ly - y_max);
      cairo_line_to(cr, x_min, ly - y_max);

      cairo_set_source_rgb(cr, r, g, b);
      cairo_fill(cr);
      cairo_set_source_rgb(cr, 0, 0, 0);
      cairo_stroke(cr);
    }
  }
  write_polygons_to_cairo_surface(cr,
                                  false,
                                  false,
                                  false,
                                  inset_state);

  if (draw_bar)
  {
    write_density_bar_to_cairo_surface(dens_min,
                                       dens_mean,
                                       dens_max,
                                       cr,
                                       bbox_bar,
                                       ly);
  }

  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

void InsetState::write_intersections_image(unsigned int res,
                                           const bool image_format_ps)
{
  std::string filename =
      inset_name() +
      "_intersections_" +
      std::to_string(n_finished_integrations());

  // Update extension
  image_format_ps ? filename += ".ps" : filename += ".svg";

  // Calculating intersections
  std::vector<Segment> intersections = intersecting_segments(res);
  cairo_surface_t *surface;

  // Create a cairo surface
  image_format_ps ? surface = cairo_ps_surface_create(filename.c_str(), lx_, ly_) : surface = cairo_svg_surface_create(filename.c_str(), lx_, ly_);
  cairo_t *cr = cairo_create(surface);

  // Write header
  if (image_format_ps)
  {

    write_ps_header(filename, surface);
  }

  write_polygons_to_cairo_surface(cr,
                                  false,
                                  false,
                                  false,
                                  this);

  cairo_set_line_width(cr, 0.0001 * std::min(lx_, ly_));

  for (auto seg : intersections)
  {

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
