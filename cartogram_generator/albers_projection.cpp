#include <fstream>
#include <iostream>

#include "cgal_typedef.h"
#include "inset_state.h"

<<<<<<< HEAD
Point albers_formula(std::vector<double> bbox, Point coords) {

  double min_lon = bbox[0];
  double min_lat = bbox[1];
  double max_lon = bbox[2];
  double max_lat = bbox[3];

  double lon = coords[0];
  double lat = coords[1];

  double radius = 1;

  // Reference Longitude and Latitude
  double lambda_0 = (min_lon + max_lon) / 2;
  double phi_0 = (min_lat + max_lat) / 2;

  // Standard Parallels
  double phi_1 = (phi_0 + max_lat) / 2;
  double phi_2 = (phi_0 + min_lat) /2;

  double n = (1/2) * (sin(phi_1) + sin(phi_2));

  double theta = n * (lon - lambda_0);
  double c = pow(cos(phi_1), 2) + (2 * n * sin(phi_1));
  double rho = (radius / n) * sqrt(c - (2 * n * sin(lat)));
  double rho_0 = (radius / n) * sqrt(c - (2 * n * sin(phi_0)));

  double new_lon = rho * sin(theta);
  double new_lat = rho_0 - (rho * cos(theta));

  Point coords_converted(new_lon, new_lat);
=======
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

Point albers_formula(CGAL::Bbox_2 bbox, Point coords) {
  // TODO
  // Convert albers_formula Python code to C++

  Point coords_converted(coords.x() + 10000, coords.y() + 10000);
>>>>>>> d6af73765ee92b5485d48d76f7644f23b54cf7bc

  return coords_converted;
}

void albers_projection(InsetState *inset_state) {
  // Get inset's bbox
  CGAL::Bbox_2 bbox = inset_bbox(inset_state);
  print_bbox(bbox);

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
<<<<<<< HEAD

  inset_state->set_geo_divs(gd_converted_vector);
}
=======
}
>>>>>>> d6af73765ee92b5485d48d76f7644f23b54cf7bc
