#include "inset_state.h"

nlohmann::json InsetState::inset_to_geojson(bool original_ext_ring_is_clockwise)
  const
{
  nlohmann::json inset_container;
  for (const auto &gd : geo_divs_) {
    nlohmann::json gd_container;
    for (const auto &pwh : gd.polygons_with_holes()) {

      // Get exterior ring of polygon with holes
      Polygon ext_ring = pwh.outer_boundary();

      // Set exterior ring to clockwise if it was originally like that
      if (original_ext_ring_is_clockwise) {
        ext_ring.reverse_orientation();
      }

      // Get exterior ring coordinates
      nlohmann::json er_container;
      for (unsigned int i = 0; i < ext_ring.size(); ++i) {
        double arr[2];
        arr[0] = ext_ring[i][0];
        arr[1] = ext_ring[i][1];
        er_container.push_back(arr);
      }

      // Repeat first point as last point as per GeoJSON standards
      er_container.push_back({ext_ring[0][0], ext_ring[0][1]});

      // Insert exterior ring into a container that stores all exterior and
      // interior rings for this polygon with holes
      nlohmann::json polygon_container;
      polygon_container.push_back(er_container);

      // Get holes of polygon with holes
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        Polygon hole = *h;

        // Set hole to counter-clockwise if it was originally like that
        if (original_ext_ring_is_clockwise) {
          hole.reverse_orientation();
        }

        // Get hole coordinates
        nlohmann::json hole_container;
        for (unsigned int i = 0; i < hole.size(); ++i) {
          double arr[2];
          arr[0] = hole[i][0];
          arr[1] = hole[i][1];
          hole_container.push_back(arr);
        }

        // Repeat first point as last point as per GeoJSON standards
        hole_container.push_back({hole[0][0], hole[0][1]});
        polygon_container.push_back(hole_container);
      }

      // Insert all polygons with holes for this GeoDiv into gd_container
      gd_container.push_back(polygon_container);
    }

    // Insert gd.id() so that we can match IDs with the coordinates when we
    // call write_to_json()
    nlohmann::json gd_id_and_coords;
    gd_id_and_coords["gd_id"] = gd.id();
    gd_id_and_coords["coordinates"] = gd_container;
    inset_container.push_back(gd_id_and_coords);
  }
  return inset_container;
}
