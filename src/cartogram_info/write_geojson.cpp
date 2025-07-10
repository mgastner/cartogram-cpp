#include "cartogram_info.hpp"
#include "constants.hpp"

// Function that returns coordinates of the end points of a "divider" line
// segment used to separate between different insets
static std::vector<double> divider_points(
  const double x1,
  const double y1,
  const double x2,
  const double y2)
{
  // Ratio between first divider point to (x2, y2) distance and
  // (x1, y1) : (x2, y2) distance
  const double ratio = divider_length + (1.0 - divider_length) / 2;

  // Calculate divider points
  const double x1d = ratio * x1 + (1.0 - ratio) * x2;
  const double x2d = ratio * x2 + (1.0 - ratio) * x1;
  const double y1d = ratio * y1 + (1.0 - ratio) * y2;
  const double y2d = ratio * y2 + (1.0 - ratio) * y1;

  // Return the two divider points (i.e., four coordinates) as a vector
  return {x1d, y1d, x2d, y2d};
}

inline std::string strip_quotes(const std::string &s)
{
  if (s.front() == '"' && s.back() == '"') {
    return s.substr(1, s.size() - 2);
  }
  return s;
}

nlohmann::json CartogramInfo::cgal_to_json(
  const bool original_geo_divs_to_geojson)
{
  nlohmann::json container = nlohmann::json::array();

  // Insert each inset into `container`
  for (const InsetState &inset_state : inset_states_) {
    std::string inset_pos = inset_state.pos();
    const nlohmann::json inset_container = inset_state.inset_to_geojson(
      original_ext_ring_is_clockwise_,
      original_geo_divs_to_geojson);

    // Insert all elements inside the inset_container (concatenate JSON arrays)
    container.insert(
      container.end(),
      inset_container.begin(),
      inset_container.end());
  }

  // Get joint bounding box for all insets.
  double bb_xmin = dbl_inf;
  double bb_ymin = dbl_inf;
  double bb_xmax = -dbl_inf;
  double bb_ymax = -dbl_inf;

  // Get the joint minimum and maximum y-coordinates for the insets "L", "C",
  // and "R". Also get the joint minimum and maximum x-coordinates for the
  // insets "T", "C", and "B". These coordinates are needed to facilitate
  // calculations for the positions of the divider lines.
  double min_xmin_tcb = dbl_inf;
  double min_ymin_lcr = dbl_inf;
  double max_xmax_tcb = -dbl_inf;
  double max_ymax_lcr = -dbl_inf;

  // Get bounding box of central inset
  Bbox inset_c_bb =
    Bbox(dbl_inf, dbl_inf, -dbl_inf, -dbl_inf);  // Initialize to invalid bbox
  for (const InsetState &inset_state : inset_states_) {
    std::string inset_pos = inset_state.pos();
    const Bbox inset_bb = inset_state.bbox(original_geo_divs_to_geojson);
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

  // Container to store divider lines
  nlohmann::json divider_container;

  // Insert divider lines between all insets
  for (const InsetState &inset_state : inset_states_) {
    std::string inset_pos = inset_state.pos();
    const Bbox inset_bb = inset_state.bbox(original_geo_divs_to_geojson);
    if (inset_pos == "T") {
      divider_container.push_back(divider_points(
        min_xmin_tcb,
        (inset_bb.ymin() + inset_c_bb.ymax()) / 2,
        max_xmax_tcb,
        (inset_bb.ymin() + inset_c_bb.ymax()) / 2));
    } else if (inset_pos == "B") {
      divider_container.push_back(divider_points(
        min_xmin_tcb,
        (inset_bb.ymax() + inset_c_bb.ymin()) / 2,
        max_xmax_tcb,
        (inset_bb.ymax() + inset_c_bb.ymin()) / 2));
    } else if (inset_pos == "L") {
      divider_container.push_back(divider_points(
        (inset_bb.xmax() + inset_c_bb.xmin()) / 2,
        max_ymax_lcr,
        (inset_bb.xmax() + inset_c_bb.xmin()) / 2,
        min_ymin_lcr));
    } else if (inset_pos == "R") {
      divider_container.push_back(divider_points(
        (inset_bb.xmin() + inset_c_bb.xmax()) / 2,
        max_ymax_lcr,
        (inset_bb.xmin() + inset_c_bb.xmax()) / 2,
        min_ymin_lcr));
    }
  }
  container.push_back(divider_container);
  return container;
}

// The format of the divider lines in the GeoJSON file is as follows:
// "dividers" : {
//   "type": "Feature",
//   "properties": { "Region": "Dividers" },
//   "geometry": {
//       "type": "MultiLineString",
//       "coordinates": [
//          [
//            [100.0, 0.0], // line1.start.x, line1.start.y
//            [101.0, 1.0] // line1.end.x, line1.end.y
//          ],
//          [
//            [102.0, 2.0],
//            [103.0, 3.0]
//          ]
//        ]
// },
static nlohmann::json add_dividers_to_geojson(const nlohmann::json &json)
{
  nlohmann::json formatted_json;
  formatted_json["type"] = "Feature";
  formatted_json["properties"]["Region"] = "Dividers";
  formatted_json["geometry"]["type"] = "MultiLineString";
  nlohmann::json coordinates = nlohmann::json::array();
  for (const auto &divider : json) {
    nlohmann::json line = nlohmann::json::array();
    line.push_back({divider[0], divider[1]});
    line.push_back({divider[2], divider[3]});
    coordinates.push_back(line);
  }
  formatted_json["geometry"]["coordinates"] = coordinates;
  return formatted_json;
}

void CartogramInfo::json_to_geojson(
  const nlohmann::json &old_json,
  nlohmann::ordered_json &new_json,
  const nlohmann::json &container)
{
  // We must match the GeoDiv IDs in the container with the IDs in the input
  // GeoJSON. For later convenience, we store the numeric indices for an ID
  // in an std::map.
  std::map<std::string, unsigned int> index_of_id_in_old_json;
  for (unsigned int index = 0; index < old_json["features"].size(); ++index) {
    const std::string id = strip_quotes(
      old_json["features"][index]["properties"].at(id_header_).dump());
    const std::pair<std::string, unsigned int> pair(id, index);
    index_of_id_in_old_json.insert(pair);
  }
  new_json["type"] = old_json["type"];
  if (n_insets() == 1) {
    new_json["bbox"] = container[(container.size() - 1)];
  } else {
    new_json["bbox"] = container[(container.size() - 2)];
    new_json["dividers"] =
      add_dividers_to_geojson(container[(container.size() - 1)]);
  }

  // Add/replace CRS to custom_crs
  new_json["crs"] = {{"type", "name"}, {"properties", {{"name", custom_crs}}}};

  new_json["properties"]["note"] =
    "Created using cartogram-cpp / go-cart.io with custom projection, not in "
    "EPSG:4326";

  // Iterate over GeoDivs and gd_ids in the container. The index
  // container.size()-2 is reserved for the bounding box, and the index
  // container.size()-1 is reserved for the divider lines. Thus, we must
  // exclude these two indices in the next loop. Hence, we only iterate over
  // n_geo_divs() elements
  for (unsigned int i = 0; i < n_geo_divs(); ++i) {
    const unsigned int index = index_of_id_in_old_json.at(
      strip_quotes(container[i].at("gd_id").dump()));
    new_json["features"][i]["type"] = "Feature";
    new_json["features"][i]["properties"] =
      old_json["features"][index]["properties"];
    new_json["features"][i]["geometry"]["type"] = "MultiPolygon";

    // Iterate over Polygon_with_holes in the GeoDiv
    for (unsigned int j = 0; j < container[i]["coordinates"].size(); ++j) {

      // Iterate over exterior ring and holes in the Polygon_with_holes
      for (unsigned int k = 0; k < container[i]["coordinates"][j].size();
           ++k) {
        new_json["features"][i]["geometry"]["coordinates"][j][k] =
          container[i]["coordinates"][j][k];
      }
    }
  }
}

void CartogramInfo::write_geojson(const std::string &suffix)
{
  std::string new_geo_file_name = map_name_ + "_" + suffix;
  std::cerr << "Writing " << new_geo_file_name << ".geojson" << std::endl;
  std::ifstream old_file(args_.geo_file_name);
  nlohmann::json old_json;
  old_file >> old_json;
  const nlohmann::json container = cgal_to_json(false);
  nlohmann::ordered_json new_json;
  json_to_geojson(old_json, new_json, container);
  if (args_.redirect_exports_to_stdout) {
    if (args_.output_equal_area_map || args_.output_shifted_insets) {
      std::cout << new_json << std::endl;
      return;
    }
    nlohmann::ordered_json new_json_original;
    nlohmann::json container_original = cgal_to_json(true);
    json_to_geojson(old_json, new_json_original, container_original);
    nlohmann::json combined_json;
    combined_json["Simplified"] = new_json;
    combined_json["Original"] = new_json_original;
    std::cout << combined_json << std::endl;
  } else {
    std::ofstream o(new_geo_file_name + ".geojson");
    o << new_json << std::endl;
  }
}