#include <iostream>
#include <fstream>
#include <regex>

#include "cgal_typedef.h"
#include "write_to_json.h"

void write_to_json(json container, std::string geo_file_name) {
  // TODO Add properties back to newJ

  std::ifstream i(geo_file_name);
  json old_j;
  i >> old_j;

  const std::string s = geo_file_name;
  std::regex rgx_geojson("[\\w_]+(?=.geojson)");
  std::regex rgx_json("[\\w_]+(?=.json)");
  std::smatch match;
  std::string gjson_fn_new = "new_geojson_simplified.json";
  if (std::regex_search(s.begin(), s.end(), match, rgx_geojson))
    gjson_fn_new = "../output_geojsons/" + match.str(0) + "_simplified.geojson";
  else if (std::regex_search(s.begin(), s.end(), match, rgx_json))
    gjson_fn_new = "../output_geojsons/" + match.str(0) + "_simplified.json";
  else
    gjson_fn_new = "new_geojson_simplified.json";

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
  
  std::ofstream o("temp.json");
  o << newJ << std::endl;

  // Replaces "aaatype" with "type" so that "type" appears at the top of the GeoJSON file
  std::ifstream in_new("temp.json");
  std::ofstream out_new(gjson_fn_new);
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
