#include "constants.hpp"
#include "inset_state.hpp"
#include "round_point.hpp"

void InsetState::adjust_for_dual_hemisphere()
{
  // Determine the maximum longitude in the western hemisphere and the minimum
  // longitude in the eastern hemisphere
  double max_lon_west = -dbl_inf;
  double min_lon_east = dbl_inf;
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto bb = pwh.bbox();
      const double xmax = bb.xmax();
      const double xmin = bb.xmin();
      max_lon_west = xmax < 0 ? std::max(xmax, max_lon_west) : max_lon_west;
      min_lon_east = xmin >= 0 ? std::min(xmin, min_lon_east) : min_lon_east;
    }
  }

  // Set transformation (translation) values to +360 for longitude
  Transformation translate(CGAL::TRANSLATION, CGAL::Vector_2<Scd>(360, 0));

  // - If min_lon_east == max_lon_west, the whole inset is contained in either
  //   only the western or only the eastern hemisphere
  // - If max_lon_west < -180.0, all polygons that are partly in the western
  //   hemisphere also are partly in the eastern hemisphere
  // - If min_lon_east > 180.0, all polygons that are partly in the eastern
  //   hemisphere also are partly in the western hemisphere
  // - If min_lon_east - max_lon_west < 180, the inset cannot fit in 1
  //   hemisphere
  if (
    max_lon_west >= -180.0 && min_lon_east <= 180.0 &&
    min_lon_east - max_lon_west >= 180) {

    // Iterate over GeoDivs
    for (auto &gd : geo_divs_) {

      // Iterate over Polygon_with_holes
      for (auto &pwh : gd.ref_to_polygons_with_holes()) {
        auto &outer_boundary = pwh.outer_boundary();

        // If pwh is in the western hemisphere
        if (pwh.bbox().xmin() < 0) {
          outer_boundary = transform(translate, outer_boundary);

          // Iterate over holes
          for (auto &h : pwh.holes()) {
            h = transform(translate, h);
          }
        }
      }
    }
  }
}

Point point_after_albers_projection(
  const Point &coords,
  double lambda_0,
  double phi_0,
  double phi_1,
  double phi_2)
{
  const double lon_in_radians = (coords.x() * pi) / 180;
  const double lat_in_radians = (coords.y() * pi) / 180;
  double x, y;
  if (abs(phi_1 + phi_2) < 1e-6) {

    // If n = 0 (i.e., phi_1 = -phi_2), the Albers projection becomes a
    // cylindrical equal-area projection with standard parallel phi_1. The
    // formula is at:
    // https://en.wikipedia.org/wiki/Cylindrical_equal-area_projection
    x = (lon_in_radians - lambda_0) * cos(phi_1);
    y = sin(lat_in_radians) / cos(phi_1);
  } else {

    // Albers projection formula:
    // https://en.wikipedia.org/wiki/Albers_projection
    const double n = 0.5 * (sin(phi_1) + sin(phi_2));
    const double c = cos(phi_1) * cos(phi_1) + 2 * n * sin(phi_1);
    const double rho_0 = sqrt(c - 2 * n * sin(phi_0)) / n;
    const double theta = n * (lon_in_radians - lambda_0);
    const double rho = sqrt(c - (2 * n * sin(lat_in_radians))) / n;
    x = rho * sin(theta);
    y = rho_0 - (rho * cos(theta));
  }
  return rounded_point({x, y}, 15);
}

void InsetState::apply_albers_projection()
{
  // Adjust the longitude coordinates if the inset spans both the eastern and
  // western hemispheres
  adjust_for_dual_hemisphere();

  // Recalculate the bbox after dual hemisphere adjustment
  const auto bb = bbox();

  // Declarations for albers_formula()
  const double min_lon = (bb.xmin() * pi) / 180;
  const double min_lat = (bb.ymin() * pi) / 180;
  const double max_lon = (bb.xmax() * pi) / 180;
  const double max_lat = (bb.ymax() * pi) / 180;

  // Reference longitude and latitude
  const double lambda_0 = 0.5 * (min_lon + max_lon);
  const double phi_0 = 0.5 * (min_lat + max_lat);

  // Standard parallels
  const double phi_1 = 0.5 * (phi_0 + max_lat);
  const double phi_2 = 0.5 * (phi_0 + min_lat);

  // Specialize/curry point_after_albers_projection() s0 that it only requires
  // one argument (Point p1).
  std::function<Point(Point)> lambda = [=](Point p1) {
    return point_after_albers_projection(p1, lambda_0, phi_0, phi_1, phi_2);
  };

  // Apply `lambda` to all points
  transform_points(lambda);
}
