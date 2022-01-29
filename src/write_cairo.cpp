#include <cairo/cairo.h>
#include <cairo/cairo-ps.h>
#include <cairo/cairo-pdf.h>
#include "constants.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include "iostream"

// TODO: IS THERE A CGAL WAY OF DETERMINING WHETHER THE LABEL'S BOUNDING
//       BOX IS COMPLETELY CONTAINED IN THE POLYGON?

// TODO: SHOULD THE CRITERION FOR PRINTING A LABEL BE THAT FITS INSIDE THE
//       POLYGON WITH HOLES? THAT CRITERION WOULD BE MORE RESTRICTIVE THAN
//       FITTING INSIDE THE EXTERIOR RING.
bool all_points_inside_exterior_ring(std::vector<Point> pts,
                                     Polygon_with_holes pwh)
{
  for (const auto &pt : pts) {
    if (pwh.outer_boundary().has_on_unbounded_side(pt)) {
      return false;
    }
  }
  return true;
}

// TODO: A SIMPLER WAY WOULD BE TO RETURN A DOUBLE. IF THE RETURN VALUE IS
//       ZERO, THEN DO NOT DRAW THE LABEL.
std::pair <bool, double> font_size(cairo_t *cr,
                                   const char *label,
                                   Point label_coordinate,
                                   GeoDiv gd)
{
  cairo_text_extents_t extents;
  double ft_size;
  for (ft_size = max_font_size; ft_size >= min_font_size; ft_size -= 0.5) {
    cairo_set_font_size(cr, ft_size);
    cairo_text_extents(cr, label, &extents);
    const auto largest_pwh = gd.largest_polygon_with_holes();

    // Bounding box of the label
    const CGAL::Bbox_2 bb(label_coordinate.x() - 0.5*extents.width,
                          label_coordinate.y() - 0.5*extents.height,
                          label_coordinate.x() + 0.5*extents.width,
                          label_coordinate.y() + 0.5*extents.height);


    // Vector of bounding-box edge points
    std::vector<Point> bb_edge_points;
    for (unsigned int i = 0; i <= 1; ++i) {
      for (unsigned int j = 0; j <= 5; ++j) {
        bb_edge_points.push_back(
          Point((j*bb.xmin() + (5 - j) * bb.xmax()) / 5,
                (i*bb.ymin() + (1 - i) * bb.ymax()))
          );
      }
    }
    if (all_points_inside_exterior_ring(bb_edge_points, largest_pwh)) {
      return std::make_pair(true, ft_size);
    }
  }
  return std::make_pair(false, 0.0);
}

void write_polygon_to_cairo_surface(cairo_t *cr,
                                    int height,
                                    bool fill_polygons,
                                    bool colors,
                                    bool plot_graticule,
                                    InsetState *inset_state)
{
  unsigned int line_width = std::min(inset_state->lx(), inset_state->ly());
  cairo_set_line_width(cr, 0.001 * line_width);
  // draw the shapes
  for (auto gd : inset_state->geo_divs()) {
    for (auto pwh : gd.polygons_with_holes()) {
      Polygon ext_ring = pwh.outer_boundary();
      // Move to starting coordinates
      cairo_move_to(cr, ext_ring[0][0], height - ext_ring[0][1]);
      // Plot each point in exterior ring
      unsigned int n = ext_ring.size();
      for (unsigned int i = 1; i < n; ++i) {
        cairo_line_to(cr, ext_ring[i][0], height - ext_ring[i][1]);
      }

      // Plot holes
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        Polygon hole = *h;
        cairo_move_to(cr, hole[0][0], height - hole[0][1]);
        unsigned int hole_size = hole.size();
        for (unsigned int i = 1; i <= hole_size; ++i) {
          cairo_line_to(cr, hole[i % hole_size][0], height - hole[i % hole_size][1]);
        }
      }

      if (colors || fill_polygons) {
        if (inset_state->is_input_target_area_missing(gd.id())) {
          // Fill path with dark-grey
          cairo_set_source_rgb(cr, 0.9375, 0.9375, 0.9375);
        } else if (colors) {
          // Get color
          Color col = inset_state->colors_at(gd.id());
          double red = col.r / 255.0;
          double green = col.g / 255.0;
          double blue = col.b / 255.0;

          // Fill path
          cairo_set_source_rgb(cr, red, green, blue);
        } else if (fill_polygons) {
          // Fill path with default color
          cairo_set_source_rgb(cr, 0.96, 0.92, 0.70);
        }
      }
      cairo_fill_preserve(cr);
      cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
      cairo_stroke(cr);
    }
  }

  //add the labels
  for (auto gd : inset_state->geo_divs()) {
    cairo_text_extents_t extents;
    //get the label
    std::string label = inset_state->labels_at(gd.id());
    const char* label_char = label.c_str();
    // go to a specific coordinate to place the label
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_select_font_face(cr, "sans-serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    Point label_coordinate = gd.point_on_surface_of_geodiv();

    // get the size of the label
    std::pair <bool, double> font_data = font_size(cr, label_char, label_coordinate, gd);
    if (font_data.first == true) {
      double font_size = font_data.second;
      cairo_set_font_size(cr, font_size);
      cairo_text_extents(cr, label_char, &extents);
      double x = label_coordinate.x() - (extents.width / 2 + extents.x_bearing);
      double y = height - label_coordinate.y() - (extents.height / 2 + extents.y_bearing);
      cairo_move_to(cr, x, y);

      cairo_show_text(cr, label_char);
    }
  }

  //plot the graticule
  if (plot_graticule) {
    boost::multi_array<XYPoint, 2> &cum_proj =
      *inset_state->ref_to_cum_proj();
    unsigned int graticule_line_spacing = 7;

    // Set line width of graticule lines
    cairo_set_line_width(cr, 0.0005 * line_width);
    cairo_set_source_rgb(cr, 0, 0, 0);

    // Vertical graticule lines
    for (unsigned int i = 0; i <= inset_state->lx(); i += graticule_line_spacing) {
      cairo_move_to(cr, cum_proj[i][0].x, cum_proj[i][0].y);
      for (unsigned int j = 1; j < inset_state->ly(); ++j) {
        cairo_line_to(cr, cum_proj[i][j].x, cum_proj[i][j].y);
      }
      cairo_stroke(cr);
    }

    // Horizontal graticule lines
    for (unsigned int j = 0; j <= inset_state->ly(); j += graticule_line_spacing) {
      cairo_move_to(cr, cum_proj[0][j].x, cum_proj[0][j].y);
      for (unsigned int i = 1; i < inset_state->lx(); ++i) {
        cairo_line_to(cr, cum_proj[i][j].x, cum_proj[i][j].y);
      }
      cairo_stroke(cr);
    }
  }
}

// Outputs a PNG file
void write_cairo_polygons_to_png(std::string fname,
                                 bool fill_polygons,
                                 bool colors,
                                 bool plot_graticule,
                                 InsetState *inset_state)
{
  const char* filename = fname.c_str();
  unsigned int width = inset_state->lx();
  unsigned int height = inset_state->ly();

  cairo_surface_t *surface;
  cairo_t *cr;

  surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
  cr = cairo_create(surface);

  write_polygon_to_cairo_surface(cr, height, fill_polygons, colors, plot_graticule, inset_state);

  cairo_surface_write_to_png(surface, filename);
  cairo_destroy(cr);
  cairo_surface_destroy(surface);
}

// Outputs a PS file
void write_cairo_polygons_to_ps(std::string fname,
                                bool fill_polygons,
                                bool colors,
                                bool plot_graticule,
                                InsetState *inset_state)
{
  const char* filename = fname.c_str();
  unsigned int width = inset_state->lx();
  unsigned int height = inset_state->ly();

  cairo_surface_t *surface;
  cairo_t *cr;

  surface = cairo_ps_surface_create(filename, width, height);
  cr = cairo_create(surface);

  //add comments
  std::string title = "%%Title: " + fname;
  cairo_ps_surface_dsc_comment(surface, title.c_str());
  cairo_ps_surface_dsc_comment(surface, "%%Creator: Michael T. Gastner et al.");
  cairo_ps_surface_dsc_comment(surface, "%%For: Humanity");
  cairo_ps_surface_dsc_comment(surface, "%%Copyright: License CC BY-NC-ND 2.0");
  cairo_ps_surface_dsc_comment(surface, "%%Magnification: 1.0000");

  write_polygon_to_cairo_surface(cr, height, fill_polygons, colors, plot_graticule, inset_state);

  cairo_show_page(cr);

  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

// Outputs both png and ps files
void write_cairo_map(std::string filename, bool plot_graticule, InsetState* inset_state) {
  std::string png_filename = filename + ".png";
  std::string ps_filename = filename + ".ps";

  //Check whether the has all GeoDivs colored
  bool has_colors = (inset_state->colors_size() == inset_state->n_geo_divs());
  write_cairo_polygons_to_png(png_filename,
                              true,
                              has_colors,
                              plot_graticule,
                              inset_state);
  write_cairo_polygons_to_ps(ps_filename,
                             true,
                             has_colors,
                             plot_graticule,
                             inset_state);
}
