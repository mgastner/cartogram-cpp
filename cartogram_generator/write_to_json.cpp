#include "cgal_typedef.h"
#include "write_to_json.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include <iostream>
#include <fstream>

json cgal_to_json(InsetState *inset_state)
{
  json container;

  // Initialize bounding box of map with bounding box of 0-th
  // Polygon_with_holes in 0-th GeoDiv
  GeoDiv gd0 = inset_state->geo_divs()[0];
  std::vector<Polygon_with_holes> pwhs = gd0.polygons_with_holes();
  CGAL::Bbox_2 bb0 = pwhs[0].bbox();
  double map_xmin = bb0.xmin();
  double map_xmax = bb0.xmax();
  double map_ymin = bb0.ymin();
  double map_ymax = bb0.ymax();

  for (auto gd : inset_state->geo_divs()) {
    json gd_container;
    for (auto pwh : gd.polygons_with_holes()) {

      // Get exterior ring of Polygon_with_holes
      Polygon ext_ring = pwh.outer_boundary();
      json polygon_container;
      json er_container;

      // Get the bounding box coordinates
      CGAL::Bbox_2 bb = pwh.bbox();
      map_xmin = (bb.xmin() < map_xmin ? bb.xmin() : map_xmin);
      map_ymin = (bb.ymin() < map_ymin ? bb.ymin() : map_ymin);
      map_xmax = (bb.xmax() > map_xmax ? bb.xmax() : map_xmax);
      map_ymax = (bb.ymax() > map_ymax ? bb.ymax() : map_ymax);

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
  container.push_back({map_xmin,map_ymin,map_xmax,map_ymax});
  return container;
}

void write_to_json(json container,
                   std::string old_geo_fn,
                   std::string new_geo_fn)
{
  // TODO Add properties back to newJ
  std::ifstream i(old_geo_fn);
  json old_j;
  i >> old_j;
  json newJ;

  // For each multipolygon in the container
  for (int i = 0; i < (int) container.size() - 1; i++) {
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
  int bbox_index = (int)container.size() - 1;
  double bbox_xmin = container[bbox_index][0];
  double bbox_ymin = container[bbox_index][1];
  double bbox_xmax = container[bbox_index][2];
  double bbox_ymax = container[bbox_index][3];
  newJ.push_back({"type", old_j["type"]});
  newJ.push_back({"bbox", {bbox_xmin, bbox_ymin, bbox_xmax, bbox_ymax}});
  std::ofstream o(new_geo_fn);
  o << newJ << std::endl;
  return;
}
