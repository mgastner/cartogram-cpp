#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <utility>

#include "cgal_typedef.h"
#include "map_state.h"

Point albers_formula(Point coords) {
  Point coords_converted = coords;

  return coords_converted;
}

void albers_projection(std::string geo_file_name, InsetState *inset_state) {
  // Get bbox from GeoJSON
  std::ifstream in_file(geo_file_name);
  nlohmann::json j;
  in_file >> j;
  std::vector<double> bbox = j["bbox"].get<std::vector<double>>();

  std::vector<GeoDiv> inset_state_converted = inset_state->geo_divs();

  for (GeoDiv gd : inset_state_converted) {
    GeoDiv gd_converted(gd.id());
    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      // Create new outer polygon
      Polygon outer_pgn_converted;

      // Iterate through outer polygon and return new coords
      Polygon outer_pgn = pgnwh.outer_boundary();
      for (Point coords_pgn : outer_pgn) {
        Point coords_pgn_converted = albers_formula(coords_pgn);
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
          Point coords_hole_converted = albers_formula(coords_hole);
          hole_converted.push_back(coords_hole_converted);
        }
        holes_v_converted.push_back(hole_converted);
      }

      Polygon_with_holes pgnwh_converted(outer_pgn_converted,
                                         holes_v_converted.begin(),
                                         holes_v_converted.end());
    }
  }

  inset_state->set_geo_divs(inset_state_converted);
}