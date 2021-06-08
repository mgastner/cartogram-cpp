#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <utility>

#include "cgal_typedef.h"
#include "map_state.h"

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

  return coords_converted;
}

void albers_projection(std::string geo_file_name, InsetState *inset_state) {
  // Get bbox from GeoJSON
  std::ifstream in_file(geo_file_name);
  nlohmann::json j;
  in_file >> j;
  std::vector<double> bbox = j["bbox"].get<std::vector<double>>();

  // Create new vector of GeoDivs
  std::vector<GeoDiv> gd_converted_vector = {};
  for (GeoDiv gd : inset_state->geo_divs()) {
    GeoDiv gd_converted(gd.id());
    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      // Create new outer polygon
      Polygon outer_pgn_converted;

      // Iterate through outer polygon and return new coords
      Polygon outer_pgn = pgnwh.outer_boundary();
      for (Point coords_pgn : outer_pgn) {
        Point coords_pgn_converted = albers_formula(bbox, coords_pgn);
        outer_pgn_converted.push_back(coords_pgn_converted);
      }

      // Create new vector of holes
      std::vector<Polygon> holes_v_converted = {};

      // Iterate through all holes
      std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
      for (Polygon hole : holes_v) {
        // Create new hole (Polygon type)
        Polygon hole_converted;

        // Iterate through each hole and return new coords for each hole
        for (Point coords_hole : hole) {
          Point coords_hole_converted = albers_formula(bbox, coords_hole);
          hole_converted.push_back(coords_hole_converted);
        }
        holes_v_converted.push_back(hole_converted);
      }

      // Create new Polygon_with_holes
      Polygon_with_holes pgnwh_converted(outer_pgn_converted,
                                         holes_v_converted.begin(),
                                         holes_v_converted.end());

      // Add new Polygon_with_holes to new GeoDiv
      gd_converted.push_back(pgnwh_converted);
    }
    // Add new GeoDiv to new vector of GeoDivs
    gd_converted_vector.push_back(gd_converted);
  }

  inset_state->set_geo_divs(gd_converted_vector);
}
