#pragma once

#include "constants.hpp"
#include "geometry.hpp"

static double choose_font_size(
  std::string_view label,
  const Point &label_pt,
  const GeoDiv &gd,
  const InsetState &)
{
  constexpr double char_factor = 0.55;  // width
  constexpr double line_factor = 0.75;  // height

  const auto &outer = gd.largest_polygon_with_holes();

  for (double fsz = max_font_size; fsz >= min_font_size; fsz -= 0.5) {
    const double w = char_factor * fsz * static_cast<double>(label.size());
    const double h = line_factor * fsz;

    CGAL::Bbox_2 bb(
      label_pt.x() - 0.5 * w,
      label_pt.y() - 0.5 * h,
      label_pt.x() + 0.5 * w,
      label_pt.y() + 0.5 * h);

    std::vector<Point> edge_samples;
    for (unsigned i = 0; i <= 1; ++i)
      for (unsigned j = 0; j <= 5; ++j)
        edge_samples.emplace_back(
          (j * bb.xmin() + (5 - j) * bb.xmax()) / 5,
          (i * bb.ymin() + (1 - i) * bb.ymax()));

    if (all_points_inside_exterior_ring(edge_samples, outer))
      return fsz;
  }
  return 0.0;
}
