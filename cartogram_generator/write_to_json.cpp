#include "constants.h"
#include "write_to_json.h"
#include <fstream>

std::vector<double> divider_points(double x1, double y1, double x2, double y2)
{

  // Ratio between first divider point to (x2, y2) distance and
  // (x1, y1) : (x2, y2) distance
  double ratio = divider_length + (1.0 - divider_length) / 2;

  //Calculate divider points
  double x1d = ratio * x1 + (1.0 - ratio) * x2;
  double x2d = ratio * x2 + (1.0 - ratio) * x1;
  double y1d = ratio * y1 + (1.0 - ratio) * y2;
  double y2d = ratio * y2 + (1.0 - ratio) * y1;

  // Return the two divider points (i.e., four coordinates) as a vector
  return {x1d, y1d, x2d, y2d};
}

nlohmann::json cgal_to_json(CartogramInfo *cart_info)
{
  nlohmann::json container;

  // Insert GeoDiv coordinates into the container
  for (auto &[inset_pos, inset_state] : *cart_info->ref_to_inset_states()) {
    for (auto gd : inset_state.geo_divs()) {
      nlohmann::json gd_container;
      for (auto pwh : gd.polygons_with_holes()) {

        // Get exterior ring of polygon with holes
        Polygon ext_ring = pwh.outer_boundary();

        // Set exterior ring to clockwise if it was originally like that
        if (cart_info->original_ext_ring_is_clockwise()) {
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
        for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
          Polygon hole = *hci;

          // Set hole to counter-clockwise if it was originally like that
          if (cart_info->original_ext_ring_is_clockwise()) {
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
      container.push_back(gd_id_and_coords);
    }
  }

  // Get joint bounding box for all insets.
  double bbox_xmin = dbl_inf;
  double bbox_ymin = dbl_inf;
  double bbox_xmax = -dbl_inf;
  double bbox_ymax = -dbl_inf;

  // Get maximum ymax and minimum ymin for "L", "C", "R" &
  // maximum xmax and minimum xmin for "T", "C", "B" insets
  // to facilitate divider line calculations
  double min_xmin_tcb = dbl_inf;
  double min_ymin_lcr = dbl_inf;
  double max_xmax_tcb = -dbl_inf;
  double max_ymax_lcr = -dbl_inf;

  // Get central inset bbox for later use on divider lines
  CGAL::Bbox_2 inset_c_bbox;
  for (auto &[inset_pos, inset_state] : *cart_info->ref_to_inset_states()) {
    CGAL::Bbox_2 inset_bbox = inset_state.bbox();
    bbox_xmin = std::min(bbox_xmin, inset_bbox.xmin());
    bbox_ymin = std::min(bbox_ymin, inset_bbox.ymin());
    bbox_xmax = std::max(bbox_xmax, inset_bbox.xmax());
    bbox_ymax = std::max(bbox_ymax, inset_bbox.ymax());
    if (inset_pos == "T" || inset_pos == "C" || inset_pos == "B") {
      max_xmax_tcb = std::max(max_xmax_tcb, inset_bbox.xmax());
      min_xmin_tcb = std::min(min_xmin_tcb, inset_bbox.xmin());
    }
    if (inset_pos == "L" || inset_pos == "C" || inset_pos == "R") {
      min_ymin_lcr = std::min(min_ymin_lcr, inset_bbox.ymin());
      max_ymax_lcr = std::max(max_ymax_lcr, inset_bbox.ymax());
    }
    if (inset_pos == "C") {
      inset_c_bbox = inset_bbox;
    }
  }

  // Insert joint bounding box into the container as a vector with four
  // numbers
  container.push_back({bbox_xmin, bbox_ymin, bbox_xmax, bbox_ymax});

  // Container to store divider lines for go-cart.io
  nlohmann::json divider_container;

  // Insert divider lines between all insets
  for (auto &[inset_pos, inset_state] : *cart_info->ref_to_inset_states()) {
    CGAL::Bbox_2 inset_bbox = inset_state.bbox();
    if (inset_pos == "R") {
      divider_container.push_back(divider_points((inset_bbox.xmin()
                                                 + inset_c_bbox.xmax()) / 2,
                                                 max_ymax_lcr,
                                                 (inset_bbox.xmin()
                                                 + inset_c_bbox.xmax()) / 2,
                                                 min_ymin_lcr));
    } else if (inset_pos == "L") {
      divider_container.push_back(divider_points((inset_bbox.xmax()
                                                 + inset_c_bbox.xmin()) / 2,
                                                 max_ymax_lcr,
                                                 (inset_bbox.xmax()
                                                 + inset_c_bbox.xmin()) / 2,
                                                 min_ymin_lcr));
    } else if (inset_pos == "T") {
      divider_container.push_back(divider_points(min_xmin_tcb,
                                                 (inset_bbox.ymin()
                                                 + inset_c_bbox.ymax()) / 2,
                                                 max_xmax_tcb,
                                                 (inset_bbox.ymin()
                                                 + inset_c_bbox.ymax()) / 2));
    } else if (inset_pos == "B") {
      divider_container.push_back(divider_points(min_xmin_tcb,
                                                 (inset_bbox.ymax()
                                                 + inset_c_bbox.ymin()) / 2,
                                                 max_xmax_tcb,
                                                 (inset_bbox.ymax()
                                                 + inset_c_bbox.ymin()) / 2));
    }
  }
  container.push_back(divider_container);
  return container;
}

void write_to_json(nlohmann::json container,
                   std::string old_geo_file_name,
                   std::string new_geo_file_name,
                   std::ostream &new_geo_stream,
                   bool output_to_stdout,
                   CartogramInfo *cart_info)
{
  std::ifstream old_file(old_geo_file_name);
  nlohmann::json old_json;
  old_file >> old_json;
  nlohmann::json new_json;

  // We must match the GeoDiv IDs in the container with the IDs in the input
  // GeoJSON. For later convenience, we store the numeric indices for an ID
  // in an std::map.
  std::map<std::string, unsigned int> index_of_id_in_old_json;
  for (unsigned int index = 0; index < old_json["features"].size(); ++index) {
    std::string id =
      old_json["features"][index]["properties"][cart_info->id_header()];
    std::pair<std::string, unsigned int> pair(id, index);
    index_of_id_in_old_json.insert(pair);
  }

  // Iterate over GeoDivs and gd_ids in the container. The index
  // container.size()-2 is reserved for the bounding box, and the index
  // container.size()-1 is reserved for the divider lines. Thus, we must
  // exclude these two indices in the next loop.
  for (unsigned int i = 0; i < container.size() - 2; ++i) {
    unsigned int index = index_of_id_in_old_json.at(container[i]["gd_id"]);
    new_json["features"][i]["properties"] =
      old_json["features"][index]["properties"];

    // TODO: THE NEXT LINE CREATES A FIELD WITH VALUE null UNLESS THE
    // INPUT GEOJSON CONTAINED A FIELD CALLED id. THIS MAY BE NEEDED BY
    // cartogram_web, BUT, IN THE LONG RUN, WE WANT TO GET RID OF IT.
    new_json["features"][i]["id"] = old_json["features"][index]["id"];
    new_json["features"][i]["type"] = "Feature";
    new_json["features"][i]["geometry"]["type"] = "MultiPolygon";

    // Iterate over Polygon_with_holes in the GeoDiv
    for (unsigned int j = 0;
         j < container[i]["coordinates"].size();
         ++j) {

      // Iterate over exterior ring and holes in the Polygon_with_holes
      for (unsigned int k = 0;
           k < container[i]["coordinates"][j].size();
           ++k) {
        new_json["features"][i]["geometry"]["coordinates"][j][k] =
          container[i]["coordinates"][j][k];
      }
    }
  }
  new_json.push_back({"type", old_json["type"]});
  new_json.push_back({"bbox", container[(container.size() - 2)]});
  new_json.push_back({"divider_points", container[(container.size() - 1)]});

  if (output_to_stdout) {
    new_geo_stream << new_json << std::endl;
  } else {
    std::ofstream o(new_geo_file_name);
    o << new_json << std::endl;
  }
  return;
}
