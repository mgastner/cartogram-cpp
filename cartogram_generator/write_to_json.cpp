#include "cgal_typedef.h"
#include "write_to_json.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include <iostream>
#include <fstream>

json cgal_to_json(InsetState *inset_state)
{
  json container;

  for (auto gd : inset_state->geo_divs()) {
    json gd_container;
    for (auto pwh : gd.polygons_with_holes()) {

      // Get exterior ring of Polygon_with_holes
      Polygon ext_ring = pwh.outer_boundary();
      json polygon_container;
      json er_container;

      for (unsigned int i = 0; i < ext_ring.size(); i++) {

        // Get exterior ring coordinates
        double arr[2];
        arr[0] = ext_ring[i][0];
        arr[1] = ext_ring[i][1];
        er_container.push_back(arr);
      }

      /* Repeat first point as last point as per GeoJSON standards. */
      er_container.push_back({ext_ring[0][0], ext_ring[0][1]});

      polygon_container.push_back(er_container);

      // Get holes of Polygon_with_holes
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); hci++) {
        Polygon hole = *hci;
        json hole_container;
        for (unsigned int i = 0; i < hole.size(); i++) {

          // Get hole coordinates
          double arr[2];
          arr[0] = hole[i][0];
          arr[1] = hole[i][1];
          hole_container.push_back(arr);
        }

        /* Repeat first point as last point as per GeoJSON standards. */
        hole_container.push_back({hole[0][0], hole[0][1]});

        polygon_container.push_back(hole_container);
      }
      gd_container.push_back(polygon_container);
    }
    container.push_back(gd_container);
  }

  return container;
}

void write_to_json(json container,
                   std::string old_geo_fn,
                   std::string new_geo_fn,
                   CGAL::Bbox_2 bb)
{
  // TODO Add properties back to newJ
  std::ifstream i(old_geo_fn);
  json old_j;
  i >> old_j;
  json newJ;

  // For each multipolygon in the container
  for (int i = 0; i < (int) container.size(); i++) {
    newJ["features"][i]["properties"] = old_j["features"][i]["properties"];
    newJ["features"][i]["id"] = old_j["features"][i]["id"];
    newJ["features"][i]["type"] = "Feature";
    newJ["features"][i]["geometry"]["type"] = "MultiPolygon";

    // For each Polygon_with_holes in the multipolygon
    for (int j = 0; j < (int) container[i].size(); j++) {

      // For each polygon or hole in the pgnWH
      for (int k = 0; k < (int) container[i][j].size(); k++) {
        newJ["features"][i]["geometry"]["coordinates"][j][k] =
          container[i][j][k];
      }
    }
  }

  newJ.push_back({"type", old_j["type"]});
  newJ.push_back({"bbox", {bb.xmin(), bb.ymin(), bb.xmax(), bb.ymax()}});

  std::ofstream o(new_geo_fn);
  o << newJ << std::endl;
  return;
}

json cgal_to_json_all_insets(CartogramInfo *cart_info) {

  json container;

  for (auto &inset_state : *cart_info->ref_to_inset_states()) {
    for (auto gd : inset_state.geo_divs()) {
      json gd_container;
      for (auto pwh : gd.polygons_with_holes()) {

        // Get exterior ring of Polygon_with_holes
        Polygon ext_ring = pwh.outer_boundary();
        json polygon_container;
        json er_container;

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
          json hole_container;
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
      container.push_back(gd_container);
    }
  }
  return container;
}

void write_to_json_all_insets(json container,
                              std::string old_geo_fn,
                              std::string new_geo_fn)
{
  std::ifstream i(old_geo_fn);
  json old_j;
  i >> old_j;
  json newJ;

  // For each multipolygon in the container
  for (int i = 0; i < (int) container.size(); i++) {
    newJ["features"][i]["properties"] = old_j["features"][i]["properties"];
    newJ["features"][i]["id"] = old_j["features"][i]["id"];
    newJ["features"][i]["type"] = "Feature";
    newJ["features"][i]["geometry"]["type"] = "MultiPolygon";

    // For each Polygon_with_holes in the multipolygon
    for (int j = 0; j < (int) container[i].size(); j++) {

      // For each polygon or hole in the pgnWH
      for (int k = 0; k < (int) container[i][j].size(); k++) {
        newJ["features"][i]["geometry"]["coordinates"][j][k] =
          container[i][j][k];
      }
    }
  }

  newJ.push_back({"type", old_j["type"]});

  std::ofstream o(new_geo_fn);
  o << newJ << std::endl;

  return;
}

// void write_to_json_all_frames(json container,
//                               std::string old_geo_fn,
//                               std::string new_geo_fn,
//                               std::map <std::string, CGAL::Bbox_2> bbox)
// {
//   std::ifstream i(old_geo_fn);
//   json old_j;
//   i >> old_j;
//   json newJ;
//
//   // For each multipolygon in the container
//   for (int i = 0; i < (int) container.size(); i++) {
//     newJ["features"][i]["properties"] = old_j["features"][i]["properties"];
//     newJ["features"][i]["id"] = old_j["features"][i]["id"];
//     newJ["features"][i]["type"] = "Feature";
//     newJ["features"][i]["geometry"]["type"] = "MultiPolygon";
//
//     // For each Polygon_with_holes in the multipolygon
//     for (int j = 0; j < (int) container[i].size(); j++) {
//
//       // For each polygon or hole in the pgnWH
//       for (int k = 0; k < (int) container[i][j].size(); k++) {
//         newJ["features"][i]["geometry"]["coordinates"][j][k] =
//           container[i][j][k];
//       }
//     }
//   }
//
//   int start_frame = container.size();
//
//   #define b element.second
//   for(auto element : bbox) {
//     newJ["features"][start_frame]["type"] = "Feature";
//     newJ["features"][start_frame]["geometry"]["type"] = "Polygon";
//     newJ["features"][start_frame]["geometry"]["coordinates"] = {{{b.xmin(),b.ymin()},
//                                                                  {b.xmax(),b.ymin()},
//                                                                  {b.xmax(),b.ymax()},
//                                                                  {b.xmin(),b.ymax()},
//                                                                  {b.xmin(),b.ymin()}}};
//     start_frame++;
//   }
//
//   newJ.push_back({"type", old_j["type"]});
//
//   std::ofstream o(new_geo_fn);
//   o << newJ << std::endl;
//
//   return;
// }
