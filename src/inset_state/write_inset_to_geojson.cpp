#include "../inset_state.h"

nlohmann::json InsetState::inset_to_json()
{
  nlohmann::json container;
  for (const auto &gd : geo_divs_) {
    nlohmann::json gd_container;
    for (const auto &pwh : gd.polygons_with_holes()) {

      // Get exterior ring of Polygon_with_holes
      Polygon ext_ring = pwh.outer_boundary();
      nlohmann::json polygon_container;
      nlohmann::json er_container;
      for (unsigned int i = 0; i < ext_ring.size(); ++i) {

        // Get exterior ring coordinates
        double arr[2];
        arr[0] = ext_ring[i][0];
        arr[1] = ext_ring[i][1];
        er_container.push_back(arr);
      }
      polygon_container.push_back(er_container);

      // Get holes of Polygon_with_holes
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        Polygon hole = *h;
        nlohmann::json hole_container;
        for (unsigned int i = 0; i < hole.size(); ++i) {

          // Get hole coordinates
          double arr[2];
          arr[0] = hole[i][0];
          arr[1] = hole[i][1];
          hole_container.push_back(arr);
        }
        polygon_container.push_back(hole_container);
      }
      gd_container.push_back(polygon_container);
    }
    container.push_back(gd_container);
  }
  return container;
}
