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

double font_size(cairo_t *cr, const char *label, Point label_pt, GeoDiv gd)
{
  cairo_text_extents_t extents;
  double ft_size;
  for (ft_size = max_font_size; ft_size >= min_font_size; ft_size -= 0.5) {
    cairo_set_font_size(cr, ft_size);
    cairo_text_extents(cr, label, &extents);
    const auto largest_pwh = gd.largest_polygon_with_holes();

    // Bounding box of the label
    const CGAL::Bbox_2 bb(label_pt.x() - 0.5*extents.width,
                          label_pt.y() - 0.5*extents.height,
                          label_pt.x() + 0.5*extents.width,
                          label_pt.y() + 0.5*extents.height);

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
      return ft_size;
    }
  }
  return 0.0;
}

void write_polygon_to_cairo_surface(cairo_t *cr,
                                    bool fill_polygons,
                                    bool colors,
                                    bool plot_graticule,
                                    InsetState *inset_state)
{
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  cairo_set_line_width(cr, 1e-3 * std::min(lx, ly));

  // Draw the shapes
  for (const auto &gd : inset_state->geo_divs()) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto ext_ring = pwh.outer_boundary();

      // Move to starting coordinates
      cairo_move_to(cr, ext_ring[0].x(), ly - ext_ring[0].y());

      // Plot each point in exterior ring
      for (unsigned int i = 1; i < ext_ring.size(); ++i) {
        cairo_line_to(cr, ext_ring[i].x(), ly - ext_ring[i].y());
      }

      // Plot holes
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        const auto hole = *h;
        cairo_move_to(cr, hole[0].x(), ly - hole[0].y());
        const unsigned int hole_size = hole.size();
        for (unsigned int i = 1; i <= hole_size; ++i) {
          cairo_line_to(cr,
                        hole[i % hole_size].x(),
                        ly - hole[i % hole_size].y());
        }
      }
      if (colors || fill_polygons) {
        if (inset_state->is_input_target_area_missing(gd.id())) {

          // Fill path with dark gray
          cairo_set_source_rgb(cr, 0.9375, 0.9375, 0.9375);
        } else if (colors) {

          // Get color
          const Color col = inset_state->color_at(gd.id());

          // Fill path
          cairo_set_source_rgb(cr,
                               col.r / 255.0,
                               col.g / 255.0,
                               col.b / 255.0);
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

  // Add labels
  for (const auto &gd : inset_state->geo_divs()) {
    const std::string label = inset_state->label_at(gd.id());
    const char* label_char = label.c_str();

    // Go to a specific coordinate to place the label
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_select_font_face(cr,
                           "sans-serif",
                           CAIRO_FONT_SLANT_NORMAL,
                           CAIRO_FONT_WEIGHT_NORMAL);
    const auto label_pt = gd.point_on_surface_of_geodiv();

    // Get size of label
    const double ft_size = font_size(cr, label_char, label_pt, gd);
    cairo_text_extents_t extents;
    if (ft_size > 0.0) {
      cairo_set_font_size(cr, ft_size);
      cairo_text_extents(cr, label_char, &extents);
      const double x =
        label_pt.x() - (extents.width / 2 + extents.x_bearing);
      const double y =
        ly - label_pt.y() - (extents.height / 2 + extents.y_bearing);
      cairo_move_to(cr, x, y);
      cairo_show_text(cr, label_char);
    }
  }

  // Plot the graticule
  if (plot_graticule) {
    const boost::multi_array<XYPoint, 2> &cum_proj =
      *inset_state->ref_to_cum_proj();
    const unsigned int graticule_line_spacing = 7;

    // Set line width of graticule lines
    cairo_set_line_width(cr, 5e-4 * std::min(lx, ly));
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);

    // Vertical graticule lines
    for (unsigned int i = 0; i <= lx; i += graticule_line_spacing) {
      cairo_move_to(cr, cum_proj[i][0].x, ly - cum_proj[i][0].y);
      for (unsigned int j = 1; j < ly; ++j) {
        cairo_line_to(cr, cum_proj[i][j].x, ly - cum_proj[i][j].y);
      }
      cairo_stroke(cr);
    }

    // Horizontal graticule lines
    for (unsigned int j = 0; j <= ly; j += graticule_line_spacing) {
      cairo_move_to(cr, cum_proj[0][j].x, ly - cum_proj[0][j].y);
      for (unsigned int i = 1; i < lx; ++i) {
        cairo_line_to(cr, cum_proj[i][j].x, ly - cum_proj[i][j].y);
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
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  cairo_surface_t *surface;
  cairo_t *cr;
  surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, lx, ly);
  cr = cairo_create(surface);
  write_polygon_to_cairo_surface(cr,
                                 fill_polygons,
                                 colors,
                                 plot_graticule,
                                 inset_state);
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
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  cairo_surface_t *surface;
  cairo_t *cr;
  surface = cairo_ps_surface_create(filename, lx, ly);
  cr = cairo_create(surface);

  // Add comments
  std::string title = "%%Title: " + fname;
  cairo_ps_surface_dsc_comment(surface, title.c_str());
  cairo_ps_surface_dsc_comment(surface,
                               "%%Creator: Michael T. Gastner et al.");
  cairo_ps_surface_dsc_comment(surface,
                               "%%For: Humanity");
  cairo_ps_surface_dsc_comment(surface,
                               "%%Copyright: License CC BY-NC-ND 2.0");
  cairo_ps_surface_dsc_comment(surface,
                               "%%Magnification: 1.0000");
  write_polygon_to_cairo_surface(cr,
                                 fill_polygons,
                                 colors,
                                 plot_graticule,
                                 inset_state);
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

// TODO: DO WE NEED THIS FUNCTION?
// Outputs both png and ps files
void write_cairo_map(std::string file_name,
                     bool plot_graticule,
                     InsetState* inset_state)
{
  const std::string png_name = file_name + ".png";
  const std::string ps_name = file_name + ".ps";

  //Check whether the has all GeoDivs colored
  const bool has_colors =
    (inset_state->colors_size() == inset_state->n_geo_divs());
  write_cairo_polygons_to_png(png_name,
                              true,
                              has_colors,
                              plot_graticule,
                              inset_state);
  write_cairo_polygons_to_ps(ps_name,
                             true,
                             has_colors,
                             plot_graticule,
                             inset_state);
}
