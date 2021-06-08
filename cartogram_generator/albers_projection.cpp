#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <utility>

#include "cgal_typedef.h"
#include "inset_state.h"

Point albers_formula(std::vector<double> bbox, Point coords) {
  // TODO
  // Convert albers_formula Python code to C++

  Point coords_converted(coords.x() + 10000, coords.y() + 10000);

  return coords_converted;
}

void albers_projection(std::string geo_file_name, InsetState *inset_state) {
  // Get bbox from GeoJSON
  std::ifstream in_file(geo_file_name);
  nlohmann::json j;
  in_file >> j;
  std::vector<double> bbox = j["bbox"].get<std::vector<double>>();

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