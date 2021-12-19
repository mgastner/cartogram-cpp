#include "cgal_typedef.h"
#include "constants.h"
#include "inset_state.h"
#include <math.h>
#include <fstream>
#include <iostream>

void print_albers_bbox(Bbox bb)
{
  std::cerr << "Bounding box of Albers-transformed map:\n";
  std::cerr << "\tx_min = "
            << bb.xmin()
            << ", y_min = "
            << bb.ymin()
            << ", x_max = "
            << bb.xmax()
            << ", y_max = "
            << bb.ymax()
            << std::endl;
  return;
}

void adjust_for_dual_hemisphere(InsetState *inset_state)
{
  // Determine the maximum longitude in the western hemisphere and the minimum
  // longitude in the eastern hemisphere
  double max_lon_west = -dbl_inf;
  double min_lon_east = dbl_inf;
  for (GeoDiv gd : (*inset_state).geo_divs()) {
    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      Bbox pgnwh_bb = pgnwh.bbox();
      double pgnwh_bb_xmax = pgnwh_bb.xmax();
      double pgnwh_bb_xmin = pgnwh_bb.xmin();
      max_lon_west = pgnwh_bb_xmax < 0
                     ? std::max(pgnwh_bb_xmax, max_lon_west)
                     : max_lon_west;
      min_lon_east = pgnwh_bb_xmin >= 0
                     ? std::min(pgnwh_bb_xmin, min_lon_east)
                     : min_lon_east;
    }
  }

  // Set transformation (translation) values to +360 for longitude
  Transformation translate(CGAL::TRANSLATION, CGAL::Vector_2<Epick>(360, 0));

  // - If min_lon_east == max_lon_west, the whole inset is contained in either
  //   only the western or only the eastern hemisphere
  // - If max_lon_west < -180.0, all polygons that are partly in the western
  //   hemisphere also are partly in the eastern hemisphere
  // - If min_lon_east > 180.0, all polygons that are partly in the eastern
  //   hemisphere also are partly in the western hemisphere
  // - If min_lon_east - max_lon_west < 180, the inset cannot fit in 1
  //   hemisphere
  if (max_lon_west >= -180.0 &&
      min_lon_east <= 180.0 &&
      min_lon_east - max_lon_west >= 180) {

    // Iterate through GeoDivs
    for (auto &gd : *inset_state->ref_to_geo_divs()) {

      // Iterate through Polygon_with_holes
      for (auto &pgnwh : *gd.ref_to_polygons_with_holes()) {
        Polygon *outer_boundary = &pgnwh.outer_boundary();

        // If the pgnwh is in the western hemisphere
        if (pgnwh.bbox().xmin() < 0) {
          *outer_boundary = transform(translate, *outer_boundary);

          // Iterate through holes
          for (auto hole_it = pgnwh.holes_begin();
               hole_it != pgnwh.holes_end();
               ++hole_it) {
            *hole_it = transform(translate, *hole_it);
          }
        }
      }
    }
  }
  return;
}

Point projected_albers_coordinates(Point coords,
                                   double lambda_0,
                                   double phi_0,
                                   double phi_1,
                                   double phi_2)
{
  double lon_in_radians = (coords.x() * pi) / 180;
  double lat_in_radians = (coords.y() * pi) / 180;
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
    double n = 0.5 * (sin(phi_1) + sin(phi_2));
    double c = cos(phi_1)*cos(phi_1) + 2*n*sin(phi_1);
    double rho_0 = sqrt(c - 2*n*sin(phi_0)) / n;
    double theta = n * (lon_in_radians - lambda_0);
    double rho = sqrt(c - (2 * n * sin(lat_in_radians))) / n;
    x = rho * sin(theta);
    y = rho_0 - (rho * cos(theta));
  }
  Point coords_converted(x, y);
  return coords_converted;
}

void transform_to_albers_projection(InsetState *inset_state)
{
  // Adjust the longitude coordinates if the inset spans both the eastern and
  // western hemispheres
  adjust_for_dual_hemisphere(inset_state);

  // Recalculate the bbox after dual hemisphere adjustment
  Bbox bb = inset_state->bbox();

  // Declarations for albers_formula()
  double min_lon = (bb.xmin() * pi) / 180;
  double min_lat = (bb.ymin() * pi) / 180;
  double max_lon = (bb.xmax() * pi) / 180;
  double max_lat = (bb.ymax() * pi) / 180;

  // Reference Longitude and Latitude
  double lambda_0 = 0.5 * (min_lon + max_lon);
  double phi_0 = 0.5 * (min_lat + max_lat);

  // Standard parallels
  double phi_1 = 0.5 * (phi_0 + max_lat);
  double phi_2 = 0.5 * (phi_0 + min_lat);

  // Iterate through GeoDivs
  for (GeoDiv &gd : *(inset_state->ref_to_geo_divs())) {

    // Iterate through Polygon_with_holes
    for (Polygon_with_holes &pgnwh : *(gd.ref_to_polygons_with_holes())) {

      // Get outer boundary
      Polygon &outer_boundary = *(&pgnwh.outer_boundary());

      // Iterate through outer boundary's coordinates
      for (Point &coords_outer : outer_boundary) {

        // Assign outer boundary's coordinates to transformed coordinates
        coords_outer = projected_albers_coordinates(coords_outer,
                                                    lambda_0,
                                                    phi_0,
                                                    phi_1,
                                                    phi_2);
      }

      // Iterate through holes
      for (auto hole_it = pgnwh.holes_begin();
           hole_it != pgnwh.holes_end();
           ++hole_it) {
        Polygon &hole = *hole_it;

        // Iterate through hole's coordinates
        for (Point &coords_hole : hole) {

          // Assign hole's coordinates to transformed coordinates
          coords_hole = projected_albers_coordinates(coords_hole,
                                                     lambda_0,
                                                     phi_0,
                                                     phi_1,
                                                     phi_2);
        }
      }
    }
  }
  return;
}
