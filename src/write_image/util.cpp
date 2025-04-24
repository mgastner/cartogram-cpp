#include "write_image.hpp"

// ======================== Utils ========================

double font_size(
  cairo_t *cr,
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
    const CGAL::Bbox_2 bb(
      label_pt.x() - 0.5 * extents.width,
      label_pt.y() - 0.5 * extents.height,
      label_pt.x() + 0.5 * extents.width,
      label_pt.y() + 0.5 * extents.height);

    // Vector of bounding-box edge points
    std::vector<Point> bb_edge_points;
    for (unsigned int i = 0; i <= 1; ++i) {
      for (unsigned int j = 0; j <= 5; ++j) {
        bb_edge_points.emplace_back(
          (j * bb.xmin() + (5 - j) * bb.xmax()) / 5,
          (i * bb.ymin() + (1 - i) * bb.ymax()));
      }
    }
    if (all_points_inside_exterior_ring(bb_edge_points, largest_pwh)) {
      return fsize;
    }
  }
  return 0.0;
}

double get_font_size(
  cairo_t *cr,
  const char *label,
  const Point label_pt,
  const GeoDiv gd)
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

void right_aligned_text(
  cairo_t *cr,
  double number,
  double x,
  double y,
  double font_size)
{
  std::string text = std::to_string(static_cast<int>(round(number)));
  cairo_text_extents_t extents;
  cairo_text_extents(cr, text.c_str(), &extents);
  double move = extents.width;
  cairo_move_to(cr, x - move - font_size * 0.75, y);
  cairo_show_text(cr, text.c_str());
}