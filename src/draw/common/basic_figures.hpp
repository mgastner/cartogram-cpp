#pragma once

#include "canvas.hpp"
#include "colors.hpp"
#include "font.hpp"
#include "geo_div.hpp"
#include "inset_state.hpp"
#include "round_point.hpp"

inline void write_triangles(
  Canvas &cvs,
  const proj_qd &proj,
  const Color &clr,
  const unsigned int ly,
  bool draw_projected_points)
{
  cvs.set_stroke(clr, 1.5e-3 * ly);

  for (Delaunay::Finite_faces_iterator fit = proj.dt.finite_faces_begin();
       fit != proj.dt.finite_faces_end();
       ++fit)
  {
    Point p1 = fit->vertex(0)->point();
    Point p2 = fit->vertex(1)->point();
    Point p3 = fit->vertex(2)->point();

    if (draw_projected_points) {
      p1 = proj.triangle_transformation.at(p1);
      p2 = proj.triangle_transformation.at(p2);
      p3 = proj.triangle_transformation.at(p3);
    }

    cvs.move_to(p1.x(), p1.y());
    cvs.line_to(p2.x(), p2.y());
    cvs.line_to(p3.x(), p3.y());
    cvs.line_to(p1.x(), p1.y());
    cvs.stroke();
  }
}

inline void write_segments(
  Canvas &cvs,
  const std::vector<Segment> &segments,
  const Color &clr)
{
  cvs.set_stroke(clr, 0.15);

  for (const auto &seg : segments) {
    Point p1 = seg.source();
    Point p2 = seg.target();

    cvs.move_to(p1.x(), p1.y());
    cvs.line_to(p2.x(), p2.y());
    cvs.stroke();
  }
}

inline void write_point(
  Canvas &cvs,
  const Point &pt,
  const Color &clr,
  const double radius = 0.5)
{
  cvs.set_fill(clr);
  cvs.circle(pt.x(), pt.y(), radius, /*filled=*/true);
  cvs.clear_fill();
}

inline void write_point_set(
  Canvas &cvs,
  std::unordered_set<Point> point_set,
  const Color &clr,
  const double radius = 0.5)
{
  for (const Point &pt : point_set)
    write_point(cvs, pt, clr, radius);
}

inline void write_polygon_points(
  Canvas &cvs,
  Color clr,
  InsetState &inset_state)
{
  for (const auto &gd : inset_state.geo_divs()) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto &ext_ring = pwh.outer_boundary();

      for (const auto &pt : ext_ring)
        write_point(cvs, pt, clr);

      for (const auto &h : pwh.holes())
        for (const auto &pt : h)
          write_point(cvs, pt, clr);
    }
  }
}

inline void write_quadtree_rectangles(
  Canvas &cvs,
  const std::vector<Bbox> &quadtree_bboxes,
  Color clr)
{
  cvs.set_stroke(clr, 0.20);

  for (const auto &bbox : quadtree_bboxes) {
    double xmin = bbox.xmin();
    double ymin = bbox.ymin();
    double xmax = bbox.xmax();
    double ymax = bbox.ymax();

    cvs.rectangle(
      xmin,
      ymin,
      xmax - xmin,
      ymax - ymin,
      /*filled =*/false);
  }
}

inline void add_white_background(Canvas &cvs, const InsetState &inset_state)
{
  cvs.paint_white(inset_state.lx(), inset_state.ly());
}

inline void write_labels(Canvas &cvs, InsetState &inset_state)
{
  cvs.set_stroke(Color{0, 0, 0}, 1.0);  // black text

  for (const auto &gd : inset_state.geo_divs()) {
    std::string label = inset_state.label_at(gd.id());
    Point pt = gd.point_on_surface_of_geodiv();

    double fsize = choose_font_size(label, pt, gd, inset_state);
    if (fsize > 0.0)
      cvs.text(pt.x(), pt.y(), label, fsize);
  }
}

inline void write_cells(Canvas &cvs, const InsetState &inset_state)
{
  cvs.set_stroke(
    Color{0, 0, 0},
    5e-4 * std::min(inset_state.lx(), inset_state.ly()));

  for (unsigned int i = 0; i < inset_state.lx() - plotted_cell_length;
       i += plotted_cell_length) {
    for (unsigned int j = 0; j < inset_state.ly() - plotted_cell_length;
         j += plotted_cell_length) {
      const Polygon cell_edge_points = inset_state.grid_cell_edge_points(i, j);

      cvs.move_to(cell_edge_points[0].x(), cell_edge_points[0].y());
      for (unsigned int k = 1; k < cell_edge_points.size(); ++k)
        cvs.line_to(cell_edge_points[k].x(), cell_edge_points[k].y());
      cvs.stroke();
    }
  }
}

inline void write_vectors(
  Canvas &cvs,
  const std::unordered_map<Point, Vector> &vectors,
  unsigned int lx,
  unsigned int ly)
{
  cvs.set_stroke(Color{0, 0, 0}, 5e-3 * std::min(lx, ly));

  for (auto const &[p, v] : vectors) {
    double theta = std::atan2(v.y(), v.x());
    double r = std::hypot(v.x(), v.y());

    double x_end = p.x() + 2 * v.x();
    double y_end = p.y() + 2 * v.y();

    cvs.move_to(p.x(), p.y());
    cvs.line_to(x_end, y_end);

    double x1 = x_end - 0.05 * r * std::cos(theta + M_PI / 6);
    double y1 = y_end - 0.05 * r * std::sin(theta + M_PI / 6);
    double x2 = x_end - 0.05 * r * std::cos(theta - M_PI / 6);
    double y2 = y_end - 0.05 * r * std::sin(theta - M_PI / 6);

    cvs.line_to(x1, y1);
    cvs.move_to(x_end, y_end);
    cvs.line_to(x2, y2);
    cvs.stroke();
  }
}

inline void write_grid(Canvas &cvs, const InsetState &inset_state)
{
  double cell = std::sqrt(compute_per_grid_cell(inset_state));
  cvs.set_stroke(
    Color{0, 0, 0},
    5e-4 * std::min(inset_state.lx(), inset_state.ly()));

  for (double x = 0; x <= inset_state.lx(); x += cell) {
    cvs.move_to(x, 0);
    cvs.line_to(x, inset_state.ly());
    cvs.stroke();
  }
  for (double y = 0; y <= inset_state.ly(); y += cell) {
    cvs.move_to(0, y);
    cvs.line_to(inset_state.lx(), y);
    cvs.stroke();
  }
}

inline void write_polygons(
  Canvas &cvs,
  bool fill_polygons,
  bool colours,
  const InsetState &inset_state,
  double line_w = 0.0,
  Color outline = Color{0, 0, 0})
{
  if (almost_equal(line_w, 0.0))
    line_w = 1e-3 * std::min(inset_state.lx(), inset_state.ly());

  cvs.set_stroke(outline, line_w);

  for (const auto &gd : inset_state.geo_divs()) {
    Color fill = outline;
    if (colours) {
      auto c = inset_state.color_at(gd.id());
      fill = Color(c.r / 255.0, c.g / 255.0, c.b / 255.0);
    } else if (fill_polygons) {
      fill = Color{0.96, 0.92, 0.70};
    }

    for (const auto &pwh : gd.polygons_with_holes()) {
      bg::model::polygon<CanvasPoint> poly;

      poly.outer().reserve(pwh.outer_boundary().size());
      for (auto const &pt : pwh.outer_boundary())
        poly.outer().push_back({pt.x(), pt.y()});

      poly.inners().resize(pwh.holes().size());
      for (std::size_t h = 0; h < pwh.holes().size(); ++h) {
        poly.inners()[h].reserve(pwh.holes()[h].size());
        for (auto const &pt : pwh.holes()[h])
          poly.inners()[h].push_back({pt.x(), pt.y()});
      }

      cvs.set_fill(fill);
      cvs.fill_polygon(poly);
      cvs.set_stroke(outline, line_w);
      cvs.stroke_polygon_outline(poly);
    }
  }

  cvs.set_stroke(Color{0, 0, 0}, 1.0);
  for (const auto &gd : inset_state.geo_divs()) {
    auto lbl = inset_state.label_at(gd.id());
    Point pt = gd.point_on_surface_of_geodiv();
    double fs = choose_font_size(lbl, pt, gd, inset_state);
    if (fs > 0)
      cvs.text(pt.x(), pt.y(), lbl, fs);
  }
}
