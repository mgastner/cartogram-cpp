#include "inset_state.hpp"
#include "write_image.hpp"
#include <unordered_set>

// ======================== Basic Plotting ========================

void write_triangles_on_surface(
  cairo_t *cr,
  const proj_qd &proj,
  const Color &clr,
  const unsigned int ly,
  bool draw_projected_points)
{
  // Draw the triangles
  for (Delaunay::Finite_faces_iterator fit = proj.dt.finite_faces_begin();
       fit != proj.dt.finite_faces_end();
       ++fit) {
    Point p1 = fit->vertex(0)->point();
    Point p2 = fit->vertex(1)->point();
    Point p3 = fit->vertex(2)->point();

    if (draw_projected_points) {
      p1 = proj.triangle_transformation.at(p1);
      p2 = proj.triangle_transformation.at(p2);
      p3 = proj.triangle_transformation.at(p3);
    }

    // set color and line width
    cairo_set_source_rgb(cr, clr.r, clr.g, clr.b);
    cairo_set_line_width(cr, 1.5e-3 * ly);

    // Start at one p1 and draw the triangle
    cairo_move_to(cr, p1.x(), ly - p1.y());
    cairo_line_to(cr, p2.x(), ly - p2.y());
    cairo_line_to(cr, p3.x(), ly - p3.y());
    cairo_line_to(cr, p1.x(), ly - p1.y());
    cairo_stroke(cr);
  }
}

// Write segments on surface
void write_segments_on_surface(
  cairo_t *cr,
  const std::vector<Segment> &segments,
  const Color &clr,
  const unsigned int ly)
{
  // Draw the segments
  for (const auto &seg : segments) {
    Point p1 = seg.source();
    Point p2 = seg.target();

    // set width of line
    cairo_set_line_width(cr, 0.15);

    // set color
    cairo_set_source_rgb(cr, clr.r, clr.g, clr.b);

    cairo_move_to(cr, p1.x(), ly - p1.y());
    cairo_line_to(cr, p2.x(), ly - p2.y());
    cairo_stroke(cr);
  }
}

void write_point_on_surface(
  cairo_t *cr,
  const Point &pt,
  const Color &clr,
  const unsigned int ly,
  const double radius = 0.5)
{
  cairo_set_source_rgb(cr, clr.r, clr.g, clr.b);
  cairo_arc(cr, pt.x(), ly - pt.y(), radius, 0, 2 * pi);
  cairo_fill(cr);
  cairo_stroke(cr);
}

void write_point_set_on_surface(
  cairo_t *cr,
  std::unordered_set<Point> point_set,
  const Color &clr,
  const unsigned int ly,
  const double radius = 0.5)
{
  for (const Point &pt : point_set) {
    write_point_on_surface(cr, pt, clr, ly, radius);
  }
}

void write_polygon_points_on_surface(
  cairo_t *cr,
  Color clr,
  InsetState &inset_state)
{
  cairo_set_source_rgb(cr, clr.r, clr.g, clr.b);
  cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
  cairo_set_line_width(cr, 0.5);
  // Draw the shapes
  for (const auto &gd : inset_state.geo_divs()) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto &ext_ring = pwh.outer_boundary();

      // Plot each point in exterior ring
      for (const auto &pt : ext_ring) {
        write_point_on_surface(cr, pt, clr, inset_state.ly());
      }

      // Plot holes
      for (const auto &h : pwh.holes()) {
        for (const auto &pt : h) {
          write_point_on_surface(cr, pt, clr, inset_state.ly());
        }
      }
    }
  }
}

void write_quadtree_rectangles_on_surface(
  cairo_t *cr,
  const std::vector<Bbox> &quadtree_bboxes,
  const Color &clr,
  const unsigned int ly)
{
  cairo_set_source_rgb(cr, clr.r, clr.g, clr.b);
  cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
  cairo_set_line_width(cr, 0.20);
  for (auto &bbox : quadtree_bboxes) {
    auto xmin = bbox.xmin();
    auto ymin = bbox.ymin();
    auto xmax = bbox.xmax();
    auto ymax = bbox.ymax();

    // draw a rectangle with bbox values
    cairo_rectangle(cr, xmin, ly - ymin, xmax - xmin, ymin - ymax);
    cairo_stroke(cr);
  }
}

// Paints cairo surface background to white
void add_white_background(cairo_t *cr)
{
  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_paint(cr);
}

void InsetState::write_delaunay_triangles(
  const std::string &filename,
  const bool draw_projected_points)
{

  std::cerr << "Writing " << filename << ".svg" << std::endl;
  cairo_surface_t *surface =
    cairo_svg_surface_create((filename + ".svg").c_str(), lx_, ly_);
  cairo_t *cr = cairo_create(surface);
  add_white_background(cr);

  write_segments_on_surface(
    cr,
    failed_constraints_dt_projected_,
    Color{1.0, 0.0, 0.0},
    ly_);
  write_segments_on_surface(
    cr,
    failed_constraints_,
    Color{0.0, 0.0, 1.0},
    ly_);
  write_triangles_on_surface(
    cr,
    proj_qd_,
    Color{0.6, 0.6, 0.6},
    ly_,
    draw_projected_points);
  write_polygons_on_surface(
    cr,
    false,
    false,
    *this,
    0.0,
    Color{1.0, 0.0, 0.0});

  // Plot points added via densification
  if (draw_projected_points) {
    project_point_set(points_from_densification_);
    project_point_set(points_before_densification_);
  }

  write_point_set_on_surface(
    cr,
    points_before_densification_,
    Color{0.0, 0.0, 0.0},
    ly_,
    0.95);
  write_point_set_on_surface(
    cr,
    points_from_densification_,
    Color{0.149, 0.545, 0.824},
    ly_,
    0.8);

  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
  failed_constraints_dt_projected_.clear();
  failed_constraints_.clear();
}

void InsetState::write_quadtree(const std::string &filename)
{
  std::cerr << "Writing " << filename << ".svg" << std::endl;
  cairo_surface_t *surface =
    cairo_svg_surface_create((filename + ".svg").c_str(), lx_, ly_);
  cairo_t *cr = cairo_create(surface);
  write_polygons_on_surface(cr, false, false, *this);
  write_quadtree_rectangles_on_surface(
    cr,
    quadtree_bboxes_,
    Color{0.6, 0.6, 0.6},
    ly_);
  write_polygon_points_on_surface(cr, Color{0.0, 0.0, 1.0}, *this);
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

// TODO: IS THERE A CGAL WAY OF DETERMINING WHETHER THE LABEL'S BOUNDING
//       BOX IS COMPLETELY CONTAINED IN THE POLYGON?

// TODO: SHOULD THE CRITERION FOR PRINTING A LABEL BE THAT IT FITS INSIDE THE
//       POLYGON WITH HOLES? THAT CRITERION WOULD BE MORE RESTRICTIVE THAN
//       FITTING INSIDE THE EXTERIOR RING.
bool all_points_inside_exterior_ring(
  const std::vector<Point> &pts,
  const Polygon_with_holes &pwh)
{
  for (const auto &pt : pts) {
    if (pwh.outer_boundary().has_on_unbounded_side(pt)) {
      return false;
    }
  }
  return true;
}

void write_polygons_on_surface(
  cairo_t *cr,
  const bool fill_polygons,
  const bool colors,
  const InsetState &inset_state,
  const double line_width,
  const Color clr)
{
  if (line_width) {
    cairo_set_line_width(cr, line_width);
  } else {
    cairo_set_line_width(
      cr,
      1e-3 * std::min(inset_state.lx(), inset_state.ly()));
  }

  // Draw the shapes
  for (const auto &gd : inset_state.geo_divs()) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto ext_ring = pwh.outer_boundary();

      // Move to starting coordinates
      cairo_move_to(cr, ext_ring[0].x(), inset_state.ly() - ext_ring[0].y());

      // Plot each point in exterior ring
      for (unsigned int i = 1; i <= ext_ring.size(); ++i) {
        cairo_line_to(
          cr,
          ext_ring[i % ext_ring.size()].x(),
          inset_state.ly() - ext_ring[i % ext_ring.size()].y());
      }

      // Plot holes
      for (const auto &h : pwh.holes()) {
        cairo_move_to(cr, h[0].x(), inset_state.ly() - h[0].y());
        for (unsigned int i = 1; i <= h.size(); ++i) {
          cairo_line_to(
            cr,
            h[i % h.size()].x(),
            inset_state.ly() - h[i % h.size()].y());
        }
      }
      if (colors || fill_polygons) {
        // if (is_input_target_area_missing(gd.id())) {

        //   // Fill path with dark gray
        //   cairo_set_source_rgb(cr, 0.9375, 0.9375, 0.9375);
        if (colors) {
          // } else if (colors) {

          // Get color
          const auto col = inset_state.color_at(gd.id());

          // Fill path
          cairo_set_source_rgb(
            cr,
            col.r / 255.0,
            col.g / 255.0,
            col.b / 255.0);
        } else if (fill_polygons) {

          // Fill path with default color
          cairo_set_source_rgb(cr, 0.96, 0.92, 0.70);
        }
        cairo_fill_preserve(cr);
      }
      cairo_set_source_rgb(cr, clr.r, clr.g, clr.b);
      cairo_stroke(cr);
    }
  }

  // Add labels
  for (const auto &gd : inset_state.geo_divs()) {
    const auto label = inset_state.label_at(gd.id());
    const auto label_char = label.c_str();

    // Go to a specific coordinate to place the label
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_select_font_face(
      cr,
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
      const double x = label_pt.x() - (extents.width / 2 + extents.x_bearing);
      const double y = inset_state.ly() - label_pt.y() -
                       (extents.height / 2 + extents.y_bearing);
      cairo_move_to(cr, x, y);
      cairo_show_text(cr, label_char);
      cairo_stroke(cr);
    }
  }

  // Plot minimum ellipses
  // for (const auto &gd : geo_divs_) {
  //   for (const auto &ell : gd.min_ellipses()) {
  //     cairo_translate(cr, ell.center.x(), ly_ - ell.center.y());
  //     cairo_rotate(cr, -ell.theta);
  //     cairo_scale(cr, ell.semimajor, ell.semiminor);
  //     cairo_arc(
  //       cr,
  //       0, //ly_ - ell.center.x(),
  //       0, //ell.center.y(),
  //       1, //std::max(ell.semimajor, ell.semiminor),
  //       0,
  //       2 * pi);
  //     cairo_identity_matrix(cr);
  //     cairo_stroke(cr);
  //   }
  // }
}

void write_vectors_on_surface(
  cairo_t *cr,
  const std::unordered_map<Point, Vector> &vectors,
  unsigned int lx_,
  unsigned int ly_)
{

  // Set line width of vectors
  cairo_set_line_width(cr, 5e-3 * std::min(lx_, ly_));
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);

  // Plot vectors
  for (const auto &v : vectors) {
    cairo_move_to(cr, v.first.x(), ly_ - v.first.y());
    // vector angle
    const auto theta = std::atan2(v.second.y(), v.second.x());
    // vector length
    const auto r =
      std::sqrt(v.second.x() * v.second.x() + v.second.y() * v.second.y());

    // vector end point
    // const auto x = v.first.x() + r * std::cos(theta);
    // const auto y = v.first.y() + r * std::sin(theta);
    const auto x = v.first.x() + 2 * v.second.x();
    const auto y = v.first.y() + 2 * v.second.y();

    // arrow head
    const auto x1 = x - 0.05 * r * std::cos(theta + pi / 6);
    const auto y1 = y - 0.05 * r * std::sin(theta + pi / 6);
    const auto x2 = x - 0.05 * r * std::cos(theta - pi / 6);
    const auto y2 = y - 0.05 * r * std::sin(theta - pi / 6);

    cairo_line_to(cr, x, ly_ - y);
    cairo_line_to(cr, x1, ly_ - y1);
    cairo_move_to(cr, x, ly_ - y);
    cairo_line_to(cr, x2, ly_ - y2);
    cairo_stroke(cr);
    cairo_stroke(cr);
  }
}

// Plot the grid
void write_grid_on_surface(cairo_t *cr, const InsetState &inset_state)
{
  const double per_grid_cell = compute_per_grid_cell(inset_state);

  std::cerr << "Total area: " << inset_state.total_inset_area() << std::endl;
  std::cerr << "Per grid cell: " << per_grid_cell << std::endl;
  const double grid_line_spacing = sqrt(per_grid_cell);
  // const unsigned int grid_line_spacing = 7;

  // Set line width of grid lines
  cairo_set_line_width(
    cr,
    5e-4 * std::min(inset_state.lx(), inset_state.ly()));
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);

  // Vertical grid lines
  for (double i = 0; i <= inset_state.lx(); i += grid_line_spacing) {
    cairo_move_to(cr, i, 0);
    cairo_line_to(cr, i, inset_state.ly());
    cairo_stroke(cr);
  }

  // Horizontal grid lines
  for (unsigned int j = 0; j <= inset_state.ly(); j += grid_line_spacing) {
    cairo_move_to(cr, 0, j);
    cairo_line_to(cr, inset_state.lx(), j);
    cairo_stroke(cr);
  }
}

// Outputs a SVG file
void write_cairo_polygons_to_svg(
  const std::string &fname,
  const bool fill_polygons,
  const bool colors,
  const bool plot_grid,
  const bool equal_area_map,
  const std::unordered_map<Point, Vector> &vectors,
  const InsetState &inset_state)
{
  const auto filename = fname.c_str();
  cairo_surface_t *surface =
    cairo_svg_surface_create(filename, inset_state.lx(), inset_state.ly());
  cairo_t *cr = cairo_create(surface);

  write_polygons_on_surface(cr, fill_polygons, colors, inset_state);
  if (plot_grid) {
    write_grid_on_surface(cr, inset_state);
    write_legend_on_surface(cr, equal_area_map, inset_state);
  }
  write_vectors_on_surface(cr, vectors, inset_state.lx(), inset_state.ly());
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

// TODO: DO WE NEED THIS FUNCTION? WOULD IT NOT MAKE MORE SENSE TO ONLY PRINT
// FILE TYPES INDICATED BY COMMAND-LINE FLAGS?
// TODO: Name should be something like write_svg, I don't think
// it needs to be exposed that we are using cairo behind the
// scenes.
// Outputs both png and SVG files
void InsetState::write_cairo_map(
  const std::string &file_name,
  const bool plot_grid,
  const bool equal_area_map,

  // TODO: This was the past signature of the function.
  // Restoring this causes a bug. Investigate.
  // const std::unordered_map<Point, Vector> &vectors)
  const std::unordered_map<Point, Vector> vectors) const
{
  const auto svg_name = file_name + ".svg";
  std::cerr << "Writing " << svg_name << std::endl;

  // Check whether the has all GeoDivs colored
  const bool has_colors = (colors_size() == n_geo_divs());
  write_cairo_polygons_to_svg(
    svg_name,
    true,
    has_colors,
    plot_grid,
    equal_area_map,
    vectors,
    *this);
}

void write_labels_on_surface(cairo_t *cr, InsetState &inset_state)
{
  for (const auto &gd : inset_state.geo_divs()) {
    const std::string label = inset_state.label_at(gd.id());
    const char *const label_char = label.c_str();
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_select_font_face(
      cr,
      "sans-serif",
      CAIRO_FONT_SLANT_NORMAL,
      CAIRO_FONT_WEIGHT_NORMAL);
    const Point label_pt = gd.point_on_surface_of_geodiv();

    // Get size of label
    const double fsize = get_font_size(cr, label_char, label_pt, gd);
    cairo_text_extents_t extents;

    // Draw label only if appropriate size is found
    if (fsize > 0.0) {
      cairo_set_font_size(cr, fsize);
      cairo_text_extents(cr, label_char, &extents);
      const double x = label_pt.x() - (extents.width / 2 + extents.x_bearing);
      const double y = inset_state.ly() - label_pt.y() -
                       (extents.height / 2 + extents.y_bearing);
      cairo_move_to(cr, x, y);
      cairo_show_text(cr, label_char);
    }
  }
}

void write_cells_on_surface(cairo_t *cr, InsetState &inset_state)
{
  // Set line width of grid lines
  cairo_set_line_width(
    cr,
    5e-4 * std::min(inset_state.lx(), inset_state.ly()));
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);

  // Iterate over grid cells
  for (unsigned int i = 0; i < inset_state.lx() - plotted_cell_length;
       i += plotted_cell_length) {
    for (unsigned int j = 0; j < inset_state.ly() - plotted_cell_length;
         j += plotted_cell_length) {

      // Draw grid cell by connecting edge points
      const Polygon cell_edge_points = inset_state.grid_cell_edge_points(i, j);
      cairo_move_to(
        cr,
        cell_edge_points[0].x(),
        inset_state.ly() - cell_edge_points[0].y());
      for (unsigned int k = 1; k < cell_edge_points.size(); ++k) {
        cairo_line_to(
          cr,
          cell_edge_points[k].x(),
          inset_state.ly() - cell_edge_points[k].y());
      }
      cairo_stroke(cr);
    }
  }
}

void InsetState::write_intersections_image()
{
  unsigned int res = intersections_resolution;
  std::string filename = inset_name() + "_intersections_" +
                         std::to_string(n_finished_integrations()) + ".svg";

  // Calculating intersections
  std::vector<Segment> intersections = intersecting_segments(res);
  cairo_surface_t *surface;

  // Create a cairo surface
  surface = cairo_svg_surface_create(filename.c_str(), lx_, ly_);
  cairo_t *cr = cairo_create(surface);

  write_polygons_on_surface(cr, false, false, *this);

  cairo_set_line_width(cr, 0.0001 * std::min(lx_, ly_));

  for (auto seg : intersections) {
    // Move to starting coordinates
    cairo_move_to(cr, seg[0].x(), ly_ - seg[0].y());

    // Draw line
    cairo_line_to(cr, seg[1].x(), ly_ - seg[1].y());

    // line with red and stroke
    cairo_set_source_rgb(cr, 1, 0, 0);
    cairo_stroke(cr);
  }
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}