#include "cgal_typedef.h"
#include "cgal_to_json.h"

json cgal_to_json(std::vector<GeoDiv> container) {
  json json_container;
  for (GeoDiv gd : container) {
    json json_multi_polygon;
    for (Polygon_with_holes pgn_wh : gd.polygons_with_holes()) {
      json json_pgn_wh_container; // Create a list of either polygons or holes

      json json_outer_pgn;
      Polygon outer_pgn = pgn_wh.outer_boundary();
      for (Point point : outer_pgn) {
        double arr[2];
        arr[0] = point.x();
        arr[1] = point.y();
        json_outer_pgn.push_back(arr);
      }
      // Add last point (first point)
      json_outer_pgn.push_back({outer_pgn[0].x(), outer_pgn[0].y()});

      // Add outer polygon to json_pgn_wh_container
      json_pgn_wh_container.push_back(json_outer_pgn);

      std::vector<Polygon> holesV(pgn_wh.holes_begin(), pgn_wh.holes_end());
      for (Polygon hole : holesV) {
        json json_hole;
        for (Point point : hole) {
          double arr[2];
          arr[0] = point.x();
          arr[1] = point.y();
          json_hole.push_back(arr);
        }
        // Add last point (first point)
        json_hole.push_back({hole[0].x(), hole[0].y()});

        // Add holes to json_pgn_wh_container
        json_pgn_wh_container.push_back(json_hole);
      }
      json_multi_polygon.push_back(json_pgn_wh_container);
    }
    json_container.push_back(json_multi_polygon);
  }
  return json_container;
}
