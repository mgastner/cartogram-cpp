#include "inset_state.hpp"

nlohmann::json InsetState::inset_to_geojson(
  bool original_ext_ring_is_clockwise,
  bool original_geo_divs_to_geojson) const
{
  auto &geo_divs =
    original_geo_divs_to_geojson ? geo_divs_original_transformed_ : geo_divs_;
  nlohmann::json inset_container;
  for (const auto &gd : geo_divs) {
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
      for (const auto &i : ext_ring) {
        er_container.push_back(
          {CGAL::to_double(i.x()), CGAL::to_double(i.y())});
      }

      // Repeat first point as last point as per GeoJSON standards
      er_container.push_back(
        {CGAL::to_double(ext_ring[0].x()), CGAL::to_double(ext_ring[0].y())});

      // Insert exterior ring into a container that stores all exterior and
      // interior rings for this polygon with holes
      nlohmann::json polygon_container;
      polygon_container.push_back(er_container);

      // Get holes of polygon with holes
      for (const auto &h : pwh.holes()) {

        // We make a copy called `hole` of *h so that `reverse_orientation()`
        // does not change the original hole
        Polygon hole = h;

        // Set hole to counter-clockwise if it was originally like that
        if (original_ext_ring_is_clockwise) {
          hole.reverse_orientation();
        }

        // Get hole coordinates
        nlohmann::json hole_container;
        for (const auto &i : hole) {
          hole_container.push_back(
            {CGAL::to_double(i.x()), CGAL::to_double(i.y())});
        }

        // Repeat first point as last point as per GeoJSON standards
        hole_container.push_back(
          {CGAL::to_double(hole[0].x()), CGAL::to_double(hole[0].y())});
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
