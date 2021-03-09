#include <iostream>
#include <fstream>
#include "cgal_typedef.h"
#include "write_to_json.h"
#include "map_state.h"

json cgal_to_json(MapState *map_state){
  
  json container;

  for (auto gd : map_state->geo_divs()){
    // For each GeoDiv
    json gd_container;

    for (auto pwh : gd.polygons_with_holes()){
      // For each PWH

      // Get PWH exterior ring
      Polygon ext_ring = pwh.outer_boundary();
      json polygon_container;
      json er_container;

      for (unsigned int i = 0; i < ext_ring.size(); i++){
        // Get exterior ring coordinates
        double arr[2];
        arr[0] = ext_ring[i][0];
        arr[1] = ext_ring[i][1];
        er_container.push_back(arr);
      }

      /* Repeat first point as last point as per GeoJSON standards. */
      er_container.push_back({ext_ring[0][0], ext_ring[0][1]});

      polygon_container.push_back(er_container);

      // Get PWH holes
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); hci++){
        Polygon hole = *hci;
        json hole_container;

        for (unsigned int i = 0; i < hole.size(); i++){
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
                   std::string new_geo_fn) {
  
  // TODO Add properties back to newJ

  std::ifstream i(old_geo_fn);
  json old_j;
  i >> old_j;
  json newJ;
  // For each multipolygon in the container 
  for (int i = 0; i < (int) container.size(); i++){
    newJ["features"][i]["properties"] = old_j["features"][i]["properties"];
    newJ["features"][i]["id"] = old_j["features"][i]["id"];

    newJ["features"][i]["type"] = "Feature";
    newJ["features"][i]["geometry"]["type"] = "MultiPolygon";
    // For each pgnWH (container of either polygons or holes) in the multipolygon 
    for (int j = 0; j < (int) container[i].size(); j++) {
      // For each polygon or hole in the pgnWH 
      for (int k = 0; k < (int) container[i][j].size(); k++) {
        newJ["features"][i]["geometry"]["coordinates"][j][k] = container[i][j][k];
      }
    }
  }
  newJ.push_back({"aaatype", old_j["type"]});
  newJ.push_back({"bbox", old_j["bbox"]});
  
  std::ofstream o("temporary.json");
  o << newJ << std::endl;

  // Replaces "aaatype" with "type" so that "type" appears at the top of the GeoJSON file
  std::ifstream in_new("temporary.json");
  std::ofstream out_new(new_geo_fn);
  std::string line;
  while (std::getline(in_new, line)) {
    while (true) {
      size_t pos = line.find("aaatype");
      if (pos != std::string::npos) 
        line.replace(pos, 7, "type");
      else
        break;
    }
    out_new << line << std::endl;
  }
}
