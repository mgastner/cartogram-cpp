#include "cgal_typedef.h"
#include "write_to_json.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include <iostream>
#include <fstream>

nlohmann::json cgal_to_json(CartogramInfo *cart_info)
{
  nlohmann::json container;
  for (auto &inset_state : *cart_info->ref_to_inset_states()) {
    for (auto gd : inset_state.geo_divs()) {
      nlohmann::json gd_container;
      for (auto pwh : gd.polygons_with_holes()) {

        // Get exterior ring of Polygon_with_holes
        Polygon ext_ring = pwh.outer_boundary();
        nlohmann::json polygon_container;
        nlohmann::json er_container;
        for (unsigned int i = 0; i < ext_ring.size(); i++) {

          // Get exterior ring coordinates
          double arr[2];
          arr[0] = ext_ring[i][0];
          arr[1] = ext_ring[i][1];
          er_container.push_back(arr);
        }
        polygon_container.push_back(er_container);

        // Get holes of Polygon_with_holes
        for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); hci++) {
          Polygon hole = *hci;
          nlohmann::json hole_container;
          for (unsigned int i = 0; i < hole.size(); i++) {

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
      nlohmann::json gd_id_and_coords;
      gd_id_and_coords["gd_id"] = gd.id();
      gd_id_and_coords["coordinates"] = gd_container;
      container.push_back(gd_id_and_coords);
    }
  }
  return container;
}

void write_to_json(nlohmann::json container,
                   std::string old_geo_fn,
                   std::string new_geo_fn,
                   std::ostream &new_geo_stream,
                   bool output_to_stdout,
                   CartogramInfo *cart_info)
{
  std::ifstream i(old_geo_fn);
  nlohmann::json old_j;
  i >> old_j;
  nlohmann::json newJ;

  // Loop over multipolygons in the container
  for (int i = 0; i < (int) container.size(); i++) {
    for (int a = 0; a < (int) old_j["features"].size(); a++) {
      if (container[i]["gd_id"] == old_j["features"][a]["properties"][cart_info->id_header()]) {
        newJ["features"][i]["properties"] = old_j["features"][a]["properties"];
        newJ["features"][i]["id"] = old_j["features"][a]["id"];
        newJ["features"][i]["type"] = "Feature";
        newJ["features"][i]["geometry"]["type"] = "MultiPolygon";

        // loop over Polygon_with_holes in the multipolygon
        for (int j = 0; j < (int) container[i]["coordinates"].size(); j++) {

          // Loop over exterior ring and holes in the Polygon_with_holes
          for (int k = 0; k < (int) container[i]["coordinates"][j].size(); k++) {
            newJ["features"][i]["geometry"]["coordinates"][j][k] =
              container[i]["coordinates"][j][k];
          }
        }
      }
    }
  }
  newJ.push_back({"type", old_j["type"]});

  if (output_to_stdout) {
    new_geo_stream << newJ << std::endl;
  } else {
    std::ofstream o(new_geo_fn);
    o << newJ << std::endl;
  }
  return;
}
