#include <fstream>
#include <iostream>
#include <math.h>

#include "cgal_typedef.h"
#include "inset_state.h"

void print_bbox(CGAL::Bbox_2 bbox) {
  std::cout << "Bounding box:" << std::endl;
  std::cout << "lon_min: " << bbox.xmin() << std::endl;
  std::cout << "lat_min: " << bbox.ymin() << std::endl;
  std::cout << "lon_max: " << bbox.xmax() << std::endl;
  std::cout << "lat_max: " << bbox.ymax() << std::endl << std::endl;
}

CGAL::Bbox_2 inset_bbox(InsetState *inset_state) {
  double inset_xmin, inset_ymin, inset_xmax, inset_ymax;

  for (GeoDiv gd : inset_state->geo_divs()) {
    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      CGAL::Bbox_2 pgnwh_bbox = pgnwh.bbox();
      inset_xmin = !inset_xmin || pgnwh_bbox.xmin() < inset_xmin
                       ? pgnwh_bbox.xmin()
                       : inset_xmin;
      inset_ymin = !inset_ymin || pgnwh_bbox.ymin() < inset_ymin
                       ? pgnwh_bbox.ymin()
                       : inset_ymin;
      inset_xmax = !inset_xmax || pgnwh_bbox.xmax() > inset_xmax
                       ? pgnwh_bbox.xmax()
                       : inset_xmax;
      inset_ymax = !inset_ymax || pgnwh_bbox.ymax() > inset_ymax
                       ? pgnwh_bbox.ymax()
                       : inset_ymax;
    }
  }

  CGAL::Bbox_2 inset_bbox(inset_xmin, inset_ymin, inset_xmax, inset_ymax);

  return inset_bbox;
}

// Declare pi globally for use in albers_formula() and albers_projection()
double pi = M_PI;

Point albers_formula(CGAL::Bbox_2 bbox, Point coords) {
  double min_lon = (bbox.xmin() * pi) / 180;
  double min_lat = (bbox.ymin() * pi) / 180;
  double max_lon = (bbox.xmax() * pi) / 180;
  double max_lat = (bbox.ymax() * pi) / 180;

  double lon = (coords.x() * pi) / 180;
  double lat = (coords.y() * pi) / 180;

  double radius = 1;

  // Reference Longitude and Latitude
  double lambda_0 = (min_lon + max_lon) / 2;
  double phi_0 = (min_lat + max_lat) / 2;

  // Standard Parallels
  double phi_1 = (phi_0 + max_lat) / 2;
  double phi_2 = (phi_0 + min_lat) / 2;

  // sin
  double sin_phi_1 = sin(phi_1);
  double sin_phi_2 = sin(phi_2);

  // cos
  double cos_phi_1 = cos(phi_1);

  double n = ((double)1 / 2) * (sin_phi_1 + sin_phi_2);

  double theta = n * (lon - lambda_0);
  double c = pow(cos_phi_1, 2) + (2 * n * sin(phi_1));
  double rho = (radius / n) * sqrt(c - (2 * n * sin(lat)));
  double rho_0 = (radius / n) * sqrt(c - (2 * n * sin(phi_0)));

  double new_lon = rho * sin(theta);
  double new_lat = rho_0 - (rho * cos(theta));

  Point coords_converted(new_lon, new_lat);

  return coords_converted;
}

void albers_projection(InsetState *inset_state) {
  // Get inset's bbox
  CGAL::Bbox_2 bbox = inset_bbox(inset_state);
  print_bbox(bbox);

  // Marc's Edit Again

  // Iterate through GeoDivs
  for (GeoDiv &gd : *(inset_state->ref_to_geo_divs())) {
    // Iterate through Polygon_with_holes
    for (Polygon_with_holes &pgnwh : *(gd.ref_to_polygons_with_holes())) {
      // Get outer boundary
      Polygon &outer_boundary = *(&pgnwh.outer_boundary());

      // Iterate through outer boundary's coordinates
      for (Point &coords_outer : outer_boundary) {
        // Assign outer boundary's coordinates to transformed coordinates
        coords_outer = albers_formula(bbox, coords_outer);
      }

      // Iterate through holes
      for (auto hole_it = pgnwh.holes_begin(); hole_it != pgnwh.holes_end();
           hole_it++) {
        Polygon &hole = *hole_it;

        // Iterate through hole's coordinates
        for (Point &coords_hole : hole) {
          // Assign hole's coordinates to transformed coordinates
          coords_hole = albers_formula(bbox, coords_hole);
        }
      }
    }
  }
}
