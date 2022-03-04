#include "../cartogram_info.h"
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

nlohmann::json CartogramInfo::cgal_to_json()
{
  nlohmann::json container = nlohmann::json::array();

  // Insert each Inset into container
  for (const auto &[inset_pos, inset_state] : inset_states_) {
    nlohmann::json inset_container =
      inset_state.inset_to_geojson(original_ext_ring_is_clockwise_);

    // Insert all elements inside the inset_container (concatenate JSON arrays)
    container.insert(container.end(),
                     inset_container.begin(),
                     inset_container.end());
  }

  // Get joint bounding box for all insets.
  double bb_xmin = dbl_inf;
  double bb_ymin = dbl_inf;
  double bb_xmax = -dbl_inf;
  double bb_ymax = -dbl_inf;

  // Get maximum ymax and minimum ymin for "L", "C", "R" &
  // maximum xmax and minimum xmin for "T", "C", "B" insets
  // to facilitate divider line calculations
  double min_xmin_tcb = dbl_inf;
  double min_ymin_lcr = dbl_inf;
  double max_xmax_tcb = -dbl_inf;
  double max_ymax_lcr = -dbl_inf;

  // Get central inset bbox for later use on divider lines
  Bbox inset_c_bb;
  for (const auto &[inset_pos, inset_state] : inset_states_) {
    Bbox inset_bb = inset_state.bbox();
    bb_xmin = std::min(bb_xmin, inset_bb.xmin());
    bb_ymin = std::min(bb_ymin, inset_bb.ymin());
    bb_xmax = std::max(bb_xmax, inset_bb.xmax());
    bb_ymax = std::max(bb_ymax, inset_bb.ymax());
    if (inset_pos == "T" || inset_pos == "C" || inset_pos == "B") {
      max_xmax_tcb = std::max(max_xmax_tcb, inset_bb.xmax());
      min_xmin_tcb = std::min(min_xmin_tcb, inset_bb.xmin());
    }
    if (inset_pos == "L" || inset_pos == "C" || inset_pos == "R") {
      min_ymin_lcr = std::min(min_ymin_lcr, inset_bb.ymin());
      max_ymax_lcr = std::max(max_ymax_lcr, inset_bb.ymax());
    }
    if (inset_pos == "C") {
      inset_c_bb = inset_bb;
    }
  }

  // Insert joint bounding box into the container as a vector with four
  // numbers
  container.push_back({bb_xmin, bb_ymin, bb_xmax, bb_ymax});

  // Divider lines are not required if there is only one inset
  if (n_insets() == 1) {
    return container;
  }

  // Container to store divider lines for go-cart.io
  nlohmann::json divider_container;

  // Insert divider lines between all insets
  for (const auto &[inset_pos, inset_state] : inset_states_) {
    Bbox inset_bb = inset_state.bbox();
    if (inset_pos == "R") {
      divider_container.push_back(divider_points((inset_bb.xmin()
                                                 + inset_c_bb.xmax()) / 2,
                                                 max_ymax_lcr,
                                                 (inset_bb.xmin()
                                                 + inset_c_bb.xmax()) / 2,
                                                 min_ymin_lcr));
    } else if (inset_pos == "L") {
      divider_container.push_back(divider_points((inset_bb.xmax()
                                                 + inset_c_bb.xmin()) / 2,
                                                 max_ymax_lcr,
                                                 (inset_bb.xmax()
                                                 + inset_c_bb.xmin()) / 2,
                                                 min_ymin_lcr));
    } else if (inset_pos == "T") {
      divider_container.push_back(divider_points(min_xmin_tcb,
                                                 (inset_bb.ymin()
                                                 + inset_c_bb.ymax()) / 2,
                                                 max_xmax_tcb,
                                                 (inset_bb.ymin()
                                                 + inset_c_bb.ymax()) / 2));
    } else if (inset_pos == "B") {
      divider_container.push_back(divider_points(min_xmin_tcb,
                                                 (inset_bb.ymax()
                                                 + inset_c_bb.ymin()) / 2,
                                                 max_xmax_tcb,
                                                 (inset_bb.ymax()
                                                 + inset_c_bb.ymin()) / 2));
    }
  }
  return container;
}

void CartogramInfo::write_geojson(std::string old_geo_file_name,
                                  std::string new_geo_file_name,
                                  std::ostream &new_geo_stream,
                                  bool output_to_stdout)
{
  std::ifstream old_file(old_geo_file_name);
  nlohmann::json old_json;
  old_file >> old_json;
  nlohmann::json new_json;
  const nlohmann::json container = cgal_to_json();

  // We must match the GeoDiv IDs in the container with the IDs in the input
  // GeoJSON. For later convenience, we store the numeric indices for an ID
  // in an std::map.
  std::map<std::string, unsigned int> index_of_id_in_old_json;
  for (unsigned int index = 0; index < old_json["features"].size(); ++index) {
    std::string id =
      old_json["features"][index]["properties"][id_header_];
    std::pair<std::string, unsigned int> pair(id, index);
    index_of_id_in_old_json.insert(pair);
  }

  // Iterate over GeoDivs and gd_ids in the container. The index
  // container.size()-2 is reserved for the bounding box, and the index
  // container.size()-1 is reserved for the divider lines. Thus, we must
  // exclude these two indices in the next loop. Hence, we only iterate over
  // n_geo_divs() elements
  for (unsigned int i = 0; i < n_geo_divs(); ++i) {
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
  new_json.push_back({"crs", "custom"});
  new_json.push_back({"divider_points", container[(container.size() - 1)]});
  if (output_to_stdout) {
    new_geo_stream << new_json << std::endl;
  } else {
    std::ofstream o(new_geo_file_name);
    o << new_json << std::endl;
  }
  return;
}