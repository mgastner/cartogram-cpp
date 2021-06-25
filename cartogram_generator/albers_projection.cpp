#include <math.h>

#include <fstream>
#include <iostream>

#include "cgal_typedef.h"
#include "inset_state.h"

void print_bbox(CGAL::Bbox_2 bbox) {
  std::cout << "Bounding box:" << std::endl;
  std::cout << "lon_min: " << bbox.xmin() << std::endl;
  std::cout << "lat_min: " << bbox.ymin() << std::endl;
  std::cout << "lon_max: " << bbox.xmax() << std::endl;
  std::cout << "lat_max: " << bbox.ymax() << std::endl << std::endl;
}

CGAL::Bbox_2 inset_bbox(InsetState inset_state) {
  double inset_xmin, inset_ymin, inset_xmax, inset_ymax;

  for (GeoDiv gd : inset_state.geo_divs()) {
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

void connect_detached_polygons() {
  // iterate through each polygon whose boundary is affected

    // iterate through each vertex in the polygon

      // identify which coordinates need to be removed
    
      // identify which coordinates need to be reconnected

      // * can't just change the coordinates, but need to add the new coordinates to the older coordinates
}

void adjust_for_dual_hemisphere(InsetState *inset_state, double bbox_xmin,
                                double bbox_xmax) {
  // Determine the maximum longitude in the western hemisphere and the minimum
  // longitude in the eastern hemisphere
  double max_lon_west = bbox_xmin;
  double min_lon_east = bbox_xmax;
  for (GeoDiv gd : (*inset_state).geo_divs()) {
    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      double pgnwh_bbox_xmax = pgnwh.bbox().xmax();
      double pgnwh_bbox_xmin = pgnwh.bbox().xmin();
      max_lon_west = pgnwh_bbox_xmax < 0 && pgnwh_bbox_xmax > max_lon_west
                         ? pgnwh_bbox_xmax
                         : max_lon_west;
      min_lon_east = pgnwh_bbox_xmin >= 0 && pgnwh_bbox_xmin < min_lon_east
                         ? pgnwh_bbox_xmin
                         : min_lon_east;
    }
  }

  // If min_lon_east == max_lon_west, the whole inset is contained in either
  // only the western or only the eastern hemisphere

  // If min_lon_east - max_lon_west < 180, the inset cannot fit in 1
  // hemisphere

  // What should be the tolerance value here? 180?
  if (min_lon_east - max_lon_west >= 180) {
    // Iterate through GeoDivs
    for (GeoDiv &gd : *(inset_state->ref_to_geo_divs())) {
      // Iterate through Polygon_with_holes
      for (Polygon_with_holes &pgnwh : *(gd.ref_to_polygons_with_holes())) {
        // Get outer boundary
        Polygon &outer_boundary = *(&pgnwh.outer_boundary());

        // Iterate through outer boundary's coordinates
        for (Point &coords_outer : outer_boundary) {
          // Assign outer boundary's coordinates to transformed coordinates

          // Translate the min_lon_east to 0 with all remaining coordinates
          // taking reference from there
          double translated_lon = coords_outer.x() >= 0
                                      ? coords_outer.x() - min_lon_east
                                      : coords_outer.x() - min_lon_east + 360;
          coords_outer = Point(translated_lon, coords_outer.y());
        }

        // Iterate through holes
        for (auto hole_it = pgnwh.holes_begin(); hole_it != pgnwh.holes_end();
             hole_it++) {
          Polygon &hole = *hole_it;

          // Iterate through hole's coordinates
          for (Point &coords_hole : hole) {
            // Assign hole's coordinates to transformed coordinates
            double translated_lon = coords_hole.x() >= 0
                                        ? coords_hole.x() - min_lon_east
                                        : coords_hole.x() - min_lon_east + 360;
            coords_hole = Point(translated_lon, coords_hole.y());
          }
        }
      }
    }
  }
}

// Declare pi globally for use in albers_formula() and albers_projection()
double pi = M_PI;

Point albers_formula(Point coords, double n, double c, double lambda_0,
                     double radius, double rho_0) {
  double lon = (coords.x() * pi) / 180;
  double lat = (coords.y() * pi) / 180;

  double theta = n * (lon - lambda_0);
  double rho = (radius / n) * sqrt(c - (2 * n * sin(lat)));

  double new_lon = rho * sin(theta);
  double new_lat = rho_0 - (rho * cos(theta));

  Point coords_converted(new_lon, new_lat);

  return coords_converted;
}

void albers_projection(InsetState *inset_state) {
  // Get inset's bbox
  CGAL::Bbox_2 bbox = inset_bbox(*inset_state);

  // Adjust the longitude coordinates if the inset spans both the eastern and
  // western hemispheres
  adjust_for_dual_hemisphere(inset_state, bbox.xmin(), bbox.xmax());

  // Recalculate the bbox after dual hemisphere adjustment
  bbox = inset_bbox(*inset_state);

  // Declarations for albers_formula()
  double min_lon = (bbox.xmin() * pi) / 180;
  double min_lat = (bbox.ymin() * pi) / 180;
  double max_lon = (bbox.xmax() * pi) / 180;
  double max_lat = (bbox.ymax() * pi) / 180;

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
  double c = pow(cos_phi_1, 2) + (2 * n * sin(phi_1));
  double rho_0 = (radius / n) * sqrt(c - (2 * n * sin(phi_0)));

  // Iterate through GeoDivs
  for (GeoDiv &gd : *(inset_state->ref_to_geo_divs())) {
    // Iterate through Polygon_with_holes
    for (Polygon_with_holes &pgnwh : *(gd.ref_to_polygons_with_holes())) {
      // Get outer boundary
      Polygon &outer_boundary = *(&pgnwh.outer_boundary());

      // Iterate through outer boundary's coordinates
      for (Point &coords_outer : outer_boundary) {
        // Assign outer boundary's coordinates to transformed coordinates
        coords_outer =
            albers_formula(coords_outer, n, c, lambda_0, radius, rho_0);
      }

      // Iterate through holes
      for (auto hole_it = pgnwh.holes_begin(); hole_it != pgnwh.holes_end();
           hole_it++) {
        Polygon &hole = *hole_it;

        // Iterate through hole's coordinates
        for (Point &coords_hole : hole) {
          // Assign hole's coordinates to transformed coordinates
          coords_hole =
              albers_formula(coords_hole, n, c, lambda_0, radius, rho_0);
        }
      }
    }
  }
}
