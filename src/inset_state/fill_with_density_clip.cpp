#include "constants.hpp"
#include "inset_state.hpp"

Polygon clip_polygon_vertical_line(const Polygon &poly, double x, bool left)
{
  Polygon clipped;
  Point prev = poly[poly.size() - 1];

  for (unsigned int i = 0; i < poly.size(); ++i) {
    const Point curr = poly[i];
    const bool prev_inside = left ? prev.x() >= x : prev.x() <= x;
    const bool curr_inside = left ? curr.x() >= x : curr.x() <= x;

    if (curr_inside) {
      if (!prev_inside) {  // Should add intersection point
        double y = prev.y() + (curr.y() - prev.y()) * (x - prev.x()) /
                                (curr.x() - prev.x());
        clipped.push_back({x, y});
      }
      clipped.push_back(curr);
    } else if (prev_inside) {  // Implies curr is outside
      // Add intersection point
      double y = prev.y() + (curr.y() - prev.y()) * (x - prev.x()) /
                              (curr.x() - prev.x());
      clipped.push_back({x, y});
    }
    prev = curr;
  }
  return clipped;
}

Polygon clip_polygon_horizontal_line(
  const Polygon &poly,
  double y,
  bool bottom)
{
  Polygon clipped;
  Point prev = poly[poly.size() - 1];

  for (unsigned int i = 0; i < poly.size(); ++i) {
    Point curr = poly[i];
    bool prev_inside = bottom ? prev.y() >= y : prev.y() <= y;
    bool curr_inside = bottom ? curr.y() >= y : curr.y() <= y;

    if (curr_inside) {
      if (!prev_inside) {  // Should add intersection point
        double x = prev.x() + (curr.x() - prev.x()) * (y - prev.y()) /
                                (curr.y() - prev.y());
        clipped.push_back({x, y});
      }
      clipped.push_back(curr);
    } else if (prev_inside) {  // Implies curr is outside
      double x = prev.x() + (curr.x() - prev.x()) * (y - prev.y()) /
                              (curr.y() - prev.y());
      clipped.push_back({x, y});
    }
    prev = curr;
  }
  return clipped;
}

double square_poly_overlap_area(const Polygon &square, const Polygon &poly)
{
  Bbox bbox = square.bbox();
  Polygon clipped = poly;
  clipped = clip_polygon_vertical_line(clipped, bbox.xmin(), true);
  clipped = clip_polygon_vertical_line(clipped, bbox.xmax(), false);
  clipped = clip_polygon_horizontal_line(clipped, bbox.ymin(), true);
  clipped = clip_polygon_horizontal_line(clipped, bbox.ymax(), false);

  return abs(clipped.area());
}

void InsetState::fill_with_density_clip(bool plot_density)
{
  // Reset the densities
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      rho_init_(i, j) = 0;
    }
  }

  // Total area accounted for each 1x1 cell
  boost::multi_array<double, 2> area_filled(boost::extents[lx_][ly_]);

  boost::multi_array<double, 2> numer(boost::extents[lx_][ly_]);
  boost::multi_array<double, 2> denom(boost::extents[lx_][ly_]);

  // Precalculate geodiv densities for fast access later
  std::vector<double> gd_target_density;
  for (const auto &gd : geo_divs_) {
    gd_target_density.push_back(target_area_at(gd.id()) / gd.area());
  }

  for (unsigned int gd_id = 0; gd_id < geo_divs_.size(); ++gd_id) {
    for (const auto &pwh : geo_divs_[gd_id].polygons_with_holes()) {
      auto bbox = pwh.outer_boundary().bbox();

      // Iterate over all 1x1 cells that may intersect the polygon
      for (unsigned int i = std::max(0, static_cast<int>(bbox.xmin()));
           i < std::min(lx_, static_cast<unsigned int>(bbox.xmax()) + 1);
           ++i) {
        for (unsigned int j = std::max(0, static_cast<int>(bbox.ymin()));
             j < std::min(ly_, static_cast<unsigned int>(bbox.ymax()) + 1);
             ++j) {

          // Create a 1x1 cell
          Polygon cell;
          cell.push_back(Point(i, j));
          cell.push_back(Point(i + 1, j));
          cell.push_back(Point(i + 1, j + 1));
          cell.push_back(Point(i, j + 1));

          double intersect_area_pwh =
            square_poly_overlap_area(cell, pwh.outer_boundary());

          // Subtract the intersection area of the holes
          for (auto hole = pwh.holes_begin(); hole != pwh.holes_end();
               ++hole) {
            intersect_area_pwh -= square_poly_overlap_area(cell, *hole);
          }

          const double weight =
            intersect_area_pwh * area_errors_[geo_divs_[gd_id].id()];
          numer[i][j] += weight * gd_target_density[gd_id];
          denom[i][j] += weight;
          area_filled[i][j] += intersect_area_pwh;
        }
      }
    }
  }

  const double ocean_density =
    (lx_ * ly_ - total_target_area()) / (lx_ * ly_ - total_inset_area());

  const double ocean_area_error =
    abs(
      (lx_ * ly_ - total_inset_area()) / (lx_ * ly_ - total_target_area()) -
      1);

  // Assume remaining area is ocean
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      double weight = (1 - area_filled[i][j]) * ocean_area_error;
      numer[i][j] += weight * ocean_density;
      denom[i][j] += weight;
    }
  }

  // Calculate the densities
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      if (denom[i][j] > 0) {
        rho_init_(i, j) = numer[i][j] / denom[i][j];
      } else { // Ocean area error is 0
        rho_init_(i, j) = ocean_density;
      }
    }
  }

  // Determine range of densities
  auto [min_iter, max_iter] = std::minmax_element(
    rho_init_.as_1d_array(),
    rho_init_.as_1d_array() + lx_ * ly_);

  dens_min_ = *min_iter;
  dens_mean_ = ocean_density;
  dens_max_ = *max_iter;

  if (plot_density) {
    std::string file_name = inset_name_ + "_unblurred_density_" +
                            std::to_string(n_finished_integrations()) + ".svg";

    std::cerr << "Writing " << file_name << std::endl;
    write_density_image(file_name, rho_init_.as_1d_array(), false);
  }

  execute_fftw_fwd_plan();
}
