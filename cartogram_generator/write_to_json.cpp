#include "cgal_typedef.h"
#include "write_to_json.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include <iostream>
#include <fstream>

std::vector<double> divider_points(double x1, double y1, double x2, double y2)
{
  double divider_length = 0.8;

  // Ratio between first divider point to (x2, y2) distance and
  // (x1, y1) : (x2, y2) distance
  double ratio = divider_length + (1 - divider_length)/2;

  //Calculate divider points
  double x1D = ratio * x1 + (1- ratio) * x2;
  double x2D = ratio * x2 + (1- ratio) * x1;
  double y1D = ratio * y1 + (1- ratio) * y2;
  double y2D = ratio * y2 + (1- ratio) * y1;

  std::vector<double> points{x1D, y1D, x2D, y2D};

  return points;
}

nlohmann::json cgal_to_json(CartogramInfo *cart_info)
{
  nlohmann::json container, divider_container;
  double bbox_xmin = 0;
  double bbox_ymin = 0;
  double bbox_xmax = 0;
  double bbox_ymax = 0;
  for (auto &[key, inset_state] : *cart_info->ref_to_inset_states()) {
    CGAL::Bbox_2 inset_bbox = inset_state.bbox();
    bbox_xmin = std::min(bbox_xmin, inset_bbox.xmin());
    bbox_ymin = std::min(bbox_ymin, inset_bbox.ymin());
    bbox_xmax = std::max(bbox_xmax, inset_bbox.xmax());
    bbox_ymax = std::max(bbox_ymax, inset_bbox.ymax());
    if (inset_state.pos() == "R") {
      divider_container.push_back(divider_points(inset_bbox.xmin(),
                                                 inset_bbox.ymax(),
                                                 inset_bbox.xmin(),
                                                 inset_bbox.ymin()));
    } else if (inset_state.pos() == "L") {
      divider_container.push_back(divider_points(inset_bbox.xmax(),
                                                 inset_bbox.ymax(),
                                                 inset_bbox.xmax(),
                                                 inset_bbox.ymin()));
    } else if (inset_state.pos() == "T") {
      divider_container.push_back(divider_points(inset_bbox.xmin(),
                                                 inset_bbox.ymin(),
                                                 inset_bbox.xmax(),
                                                 inset_bbox.ymin()));
    } else if (inset_state.pos() == "B") {
      divider_container.push_back(divider_points(inset_bbox.xmin(),
                                                 inset_bbox.ymax(),
                                                 inset_bbox.xmax(),
                                                 inset_bbox.ymax()));
    }
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
      container.push_back(gd_container);
    }
  }
  container.push_back({bbox_xmin, bbox_ymin, bbox_xmax, bbox_ymax});
  container.push_back(divider_container);
  return container;
}

void write_to_json(nlohmann::json container,
                   std::string old_geo_fn,
                   std::string new_geo_fn,
                   std::ostream &new_geo_stream,
                   bool output_to_stdout)
{
  std::ifstream i(old_geo_fn);
  nlohmann::json old_j;
  i >> old_j;
  nlohmann::json newJ;

  // Loop over multipolygons in the container
  for (int i = 0; i < (int) container.size() -2; i++) {
    newJ["features"][i]["properties"] = old_j["features"][i]["properties"];
    newJ["features"][i]["id"] = old_j["features"][i]["id"];
    newJ["features"][i]["type"] = "Feature";
    newJ["features"][i]["geometry"]["type"] = "MultiPolygon";

    // loop over Polygon_with_holes in the multipolygon
    for (int j = 0; j < (int) container[i].size(); j++) {

      // Loop over exterior ring and holes in the Polygon_with_holes
      for (int k = 0; k < (int) container[i][j].size(); k++) {
        newJ["features"][i]["geometry"]["coordinates"][j][k] =
          container[i][j][k];
      }
    }
  }
  newJ.push_back({"type", old_j["type"]});
  newJ.push_back({"bbox", container[(container.size() - 2)]});
  newJ.push_back({"divider_points", container[(container.size() - 1)]});
  if (output_to_stdout) {
    new_geo_stream << newJ << std::endl;
  } else {
    std::ofstream o(new_geo_fn);
    o << newJ << std::endl;
  }
  return;
}
