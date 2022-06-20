#include "constants.h"
#include "inset_state.h"
#include <cairo/cairo-ps.h>
#include <cairo/cairo-pdf.h>
#include <iostream>

void write_triangles_on_cairo_surface(cairo_t *cr, Delaunay &dt, color clr)
{
  // Draw the triangles
  for (Delaunay::Finite_faces_iterator fit = dt.finite_faces_begin();
       fit != dt.finite_faces_end();
       ++fit) {
    Point p1 = fit->vertex(0)->point();
    Point p2 = fit->vertex(1)->point();
    Point p3 = fit->vertex(2)->point();
    
    // set width of line
    cairo_set_line_width(cr, 0.05);

    // set color
    cairo_set_source_rgb(cr, clr.r, clr.g, clr.b);

    cairo_move_to(cr, p1.x(), 512 - p1.y());
    cairo_line_to(cr, p2.x(), 512 - p2.y());
    cairo_line_to(cr, p3.x(), 512 - p3.y());
    cairo_line_to(cr, p1.x(), 512 - p1.y());
    cairo_stroke(cr);
  }
}

void write_ps_header(std::string filename, cairo_surface_t *surface)
{
  const std::string title = "%%Title: " + filename;
  cairo_ps_surface_dsc_comment(surface, title.c_str());
  cairo_ps_surface_dsc_comment(
    surface,
    "%%Creator: Michael T. Gastner et al.");
  cairo_ps_surface_dsc_comment(surface, "%%For: Humanity");
  cairo_ps_surface_dsc_comment(surface, "%%Copyright: License CC BY");
  cairo_ps_surface_dsc_comment(surface, "%%Magnification: 1.0000");
}

void write_point_on_cairo_surface(cairo_t *cr, Point pt, color clr)
{
  cairo_set_source_rgb(cr, clr.r, clr.g, clr.b);
  cairo_move_to(cr, pt.x(), 512 - pt.y());
  cairo_line_to(cr, pt.x() + 0.05, 512 - pt.y() - 0.05);
  cairo_close_path(cr);
  cairo_stroke(cr);;
}


void InsetState::write_polygon_points_on_cairo_surface(cairo_t *cr, color clr)
{
  cairo_set_source_rgb(cr, clr.r, clr.g, clr.b);
  cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
  cairo_set_line_width(cr, 0.5);
  // Draw the shapes
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto ext_ring = pwh.outer_boundary();

      // Plot each point in exterior ring
      for (unsigned int i = 0; i < ext_ring.size(); ++i) {
        write_point_on_cairo_surface(cr, ext_ring[i], clr);
      }

      // Plot holes
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
        Polygon hole = *hci;
        for (unsigned int i = 1; i <= hole.size(); ++i) {
          write_point_on_cairo_surface(cr, hole[i], clr);
        }
      }
    }
  }
}

void InsetState::write_quadtree(const std::string filename) {
  cairo_surface_t *surface = cairo_ps_surface_create(
    (filename + "_dt_test.ps").c_str(),
    lx_,
    ly_);
  cairo_t *cr = cairo_create(surface);
  write_ps_header((filename + "_dt_test.ps"), surface);
  write_triangles_on_cairo_surface(cr, proj_qd_.dt, color{0.0, 0.0, 0.0});
  write_polygon_points_on_cairo_surface(cr, color{1.0, 0.0, 0.0});
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

// TODO: IS THERE A CGAL WAY OF DETERMINING WHETHER THE LABEL'S BOUNDING
//       BOX IS COMPLETELY CONTAINED IN THE POLYGON?

// TODO: SHOULD THE CRITERION FOR PRINTING A LABEL BE THAT IT FITS INSIDE THE
//       POLYGON WITH HOLES? THAT CRITERION WOULD BE MORE RESTRICTIVE THAN
//       FITTING INSIDE THE EXTERIOR RING.
bool all_points_inside_exterior_ring(const std::vector<Point> pts,
                                     const Polygon_with_holes pwh)
{
  for (const auto &pt : pts) {
    if (pwh.outer_boundary().has_on_unbounded_side(pt)) {
      return false;
    }
  }
  return true;
}

double font_size(cairo_t *cr,
                 const char *label,
                 const Point label_pt,
                 const GeoDiv &gd)
{
  for (double fsize = max_font_size; fsize >= min_font_size; fsize -= 0.5) {
    cairo_set_font_size(cr, fsize);
    cairo_text_extents_t extents;
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
      return fsize;
    }
  }
  return 0.0;
}

void InsetState::write_polygons_to_cairo_surface(cairo_t *cr,
                                                 const bool fill_polygons,
                                                 const bool colors,
                                                 const bool plot_graticule)
{
  cairo_set_line_width(cr, 1e-3 * std::min(lx_, ly_));

  // Draw the shapes
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto ext_ring = pwh.outer_boundary();

      // Move to starting coordinates
      cairo_move_to(cr, ext_ring[0].x(), ly_ - ext_ring[0].y());

      // Plot each point in exterior ring
      for (unsigned int i = 1; i < ext_ring.size(); ++i) {
        cairo_line_to(cr, ext_ring[i].x(), ly_ - ext_ring[i].y());
      }

      // Plot holes
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        cairo_move_to(cr, (*h)[0].x(), ly_ - (*h)[0].y());
        const unsigned int hsize = h->size();
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
          const auto col = color_at(gd.id());

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
  for (const auto &gd : geo_divs_) {
    const auto label = label_at(gd.id());
    const auto label_char = label.c_str();

    // Go to a specific coordinate to place the label
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_select_font_face(cr,
                           "sans-serif",
                           CAIRO_FONT_SLANT_NORMAL,
                           CAIRO_FONT_WEIGHT_NORMAL);
    const auto label_pt = gd.point_on_surface_of_geodiv();

    // Get size of label
    const auto fsize = font_size(cr, label_char, label_pt, gd);
    if (fsize > 0.0) {
      cairo_set_font_size(cr, fsize);
      cairo_text_extents_t extents;
      cairo_text_extents(cr, label_char, &extents);
      const double x =
        label_pt.x() - (extents.width / 2 + extents.x_bearing);
      const double y =
        ly_ - label_pt.y() - (extents.height / 2 + extents.y_bearing);
      cairo_move_to(cr, x, y);
      cairo_show_text(cr, label_char);
    }
  }

  // Plot the graticule
  if (plot_graticule) {
    const unsigned int graticule_line_spacing = 7;

    // Set line width of graticule lines
    cairo_set_line_width(cr, 5e-4 * std::min(lx_, ly_));
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);

    // Vertical graticule lines
    for (unsigned int i = 0; i <= lx_; i += graticule_line_spacing) {
      cairo_move_to(cr, cum_proj_[i][0].x, ly_ - cum_proj_[i][0].y);
      for (unsigned int j = 1; j < ly_; ++j) {
        cairo_line_to(cr, cum_proj_[i][j].x, ly_ - cum_proj_[i][j].y);
      }
      cairo_stroke(cr);
    }

    // Horizontal graticule lines
    for (unsigned int j = 0; j <= ly_; j += graticule_line_spacing) {
      cairo_move_to(cr, cum_proj_[0][j].x, ly_ - cum_proj_[0][j].y);
      for (unsigned int i = 1; i < lx_; ++i) {
        cairo_line_to(cr, cum_proj_[i][j].x, ly_ - cum_proj_[i][j].y);
      }
      cairo_stroke(cr);
    }
  }
}

// Outputs a PNG file
void InsetState::write_cairo_polygons_to_png(const std::string fname,
                                             const bool fill_polygons,
                                             const bool colors,
                                             const bool plot_graticule)
{
  const auto filename = fname.c_str();
  cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32,
                                                        lx_,
                                                        ly_);
  cairo_t *cr = cairo_create(surface);
  write_polygons_to_cairo_surface(cr,
                                  fill_polygons,
                                  colors,
                                  plot_graticule);
  cairo_surface_write_to_png(surface, filename);
  cairo_destroy(cr);
  cairo_surface_destroy(surface);
}

// Outputs a PS file
void InsetState::write_cairo_polygons_to_ps(const std::string fname,
                                            const bool fill_polygons,
                                            const bool colors,
                                            const bool plot_graticule)
{
  const auto filename = fname.c_str();
  cairo_surface_t *surface = cairo_ps_surface_create(filename, lx_, ly_);
  cairo_t *cr = cairo_create(surface);

  // Add comments
  const std::string title = "%%Title: " + fname;
  cairo_ps_surface_dsc_comment(surface, title.c_str());
  cairo_ps_surface_dsc_comment(surface,
                               "%%Creator: Michael T. Gastner et al.");
  cairo_ps_surface_dsc_comment(surface,
                               "%%For: Humanity");
  cairo_ps_surface_dsc_comment(surface,
                               "%%Copyright: License CC BY");
  cairo_ps_surface_dsc_comment(surface,
                               "%%Magnification: 1.0000");
  write_polygons_to_cairo_surface(cr,
                                  fill_polygons,
                                  colors,
                                  plot_graticule);
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

// TODO: DO WE NEED THIS FUNCTION? WOULD IT NOT MAKE MORE SENSE TO ONLY PRINT
// FILE TYPES INDICATED BY COMMAND-LINE FLAGS?
// Outputs both png and ps files
void InsetState::write_cairo_map(const std::string file_name,
                                 const bool plot_graticule)
{
  const auto png_name = file_name + ".png";
  const auto ps_name = file_name + ".ps";

  //Check whether the has all GeoDivs colored
  const bool has_colors =
    (colors_size() == n_geo_divs());
  write_cairo_polygons_to_png(png_name,
                              true,
                              has_colors,
                              plot_graticule);
  write_cairo_polygons_to_ps(ps_name,
                             true,
                             has_colors,
                             plot_graticule);
}
