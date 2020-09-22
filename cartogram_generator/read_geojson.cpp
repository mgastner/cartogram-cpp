// TODO:
// - IMPROVE ERROR HANDLING
// - ADD SUPPORT FOR "geometry": "Polygon"
// - CHECK WHETHER COORDINATES ARE EMPTY

#include "geo_div.h"
#include <nlohmann/json.hpp>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <fstream>

using json = nlohmann::json;

void check_geojson_validity(json j)
{
  if (!j.contains(std::string{"type"})) {
    std::cerr << "ERROR: JSON does not contain a key 'type'" << std::endl;
    _Exit(4);
  }
  if (j["type"] != "FeatureCollection") {
    std::cerr << "ERROR: JSON is not a valid GeoJSON FeatureCollection"
              << std::endl;
    _Exit(5);
  }
  if (!j.contains(std::string{"features"})) {
    std::cerr << "ERROR: JSON does not contain a key 'features'" << std::endl;
    _Exit(6);
  }
  json features = j["features"];
  for (auto feature : features) {
    if (!feature.contains(std::string{"type"})) {
      std::cerr << "ERROR: JSON contains a 'Features' element without key "
                << "'type'" << std::endl;
      _Exit(7);
    }
    if (feature["type"] != "Feature") {
      std::cerr << "ERROR: JSON contains a 'Features' element whose type "
                << "is not 'Feature'" << std::endl;
      _Exit(8);
    }
    if (!feature.contains(std::string("geometry"))) {
      std::cerr << "ERROR: JSON contains a feature without key 'geometry'"
                << std::endl;
      _Exit(9);
    }
    json geometry = feature["geometry"];
    if (!geometry.contains(std::string("type"))) {
      std::cerr << "ERROR: JSON contains geometry without key 'type'"
                << std::endl;
      _Exit(10);
    }
    if (!geometry.contains(std::string("coordinates"))) {
      std::cerr << "ERROR: JSON contains geometry without key 'coordinates'"
                << std::endl;
      _Exit(11);
    }
    if (geometry["type"] != "MultiPolygon" && geometry["type"] != "Polygon") {
      std::cerr << "ERROR: JSON contains unsupported geometry "
                << geometry["type"] << std::endl;
      _Exit(12);
    }
  }

  return;
}

GeoDiv JSONToCGAL(std::string id, json json_coords) {
  GeoDiv gd(id);
  for (json json_pgn_holes_container : json_coords) {
    using namespace CGAL;
    Polygon_2<K> cgal_pgn;
    //std::vector<Polygon_2<K> > holesV;
    // Add outer polygon
    for (int j = 0; j < json_pgn_holes_container[0].size(); j++) {
      cgal_pgn.push_back(K::Point_2(json_pgn_holes_container[0][j][0],
                                    json_pgn_holes_container[0][j][1]));
    }

    set_pretty_mode(std::cout);
    std::cout << cgal_pgn << std::endl;

    //for (int i = 1; i < json_pgn_holes_container.size(); i++) {
    //Polygon_2 hole;
//     // Add hole(s)
//     for (int j = 0; j < JSONPgnWHContainer[i].size(); j++) {
//       hole.push_back(Point_2(JSONPgnWHContainer[i][j][0], JSONPgnWHContainer[i][j][1]));
//     }
//     holesV.push_back(hole);
//   }
// }
// if (holesV.empty()) {
//   PolygonWH pgnWH(CGALPgn);
//   gD.addPolygon(pgnWH);
// } else {
//   PolygonWH pgnWH(CGALPgn, holesV.begin(), holesV.end());
//   gD.addPolygon(pgnWH);
  }
  return gd;
}

void read_geojson(std::string geometry_file_name)
{
  // Open file.
  std::ifstream in_file(geometry_file_name);
  if (!in_file) {
    throw std::system_error(errno,
                            std::system_category(),
                            "failed to open " + geometry_file_name);
  }

  // Parse JSON.
  json j;
  try {
    in_file >> j;
  } catch (json::parse_error& e) {
    std::cerr << "ERROR: " << e.what() << '\n'
              << "exception id: " << e.id << '\n'
              << "byte position of error: " << e.byte << std::endl;
    _Exit(3);
  }
  check_geojson_validity(j);
  for (auto feature : j["features"]) {
    json geometry = feature["geometry"];
    if (geometry["type"] == "Polygon") {
      std::cerr << "ERROR: Sorry, no support for Polygon geometry yet"
                << std::endl;
      _Exit(13);
    } else if (geometry["type"] == "MultiPolygon") {
      JSONToCGAL("id", geometry["coordinates"]);
    }
  }
  return;
}
