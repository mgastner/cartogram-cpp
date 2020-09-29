// TODO:
// - IMPROVE ERROR HANDLING
// - ADD SUPPORT FOR "geometry": "Polygon"
// - CHECK WHETHER COORDINATES ARE EMPTY
// - ENFORCE ORIENTATION RULE FOR EXTERIOR RINGS AND INTERIOR RINGS IN read_geojson()
// - AVOID REPEATING THE STARTING POINT IN THE CGAL POLYGON.

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

    // Store exterior ring in CGAL format
    Polygon_2<K> ext_ring;
    json jphc_ext = json_pgn_holes_container[0];
    for (int j = 0; j < jphc_ext.size() - 1; j++) {
      ext_ring.push_back(K::Point_2(jphc_ext[j][0], jphc_ext[j][1]));
    }

    // CGAL considers a polygon as simple only if first vertex and last vertex
    // are different
    if (jphc_ext[0][0] != jphc_ext[jphc_ext.size() - 1][0] ||
        jphc_ext[0][1] != jphc_ext[jphc_ext.size() - 1][1]) {
      ext_ring.push_back(K::Point_2(jphc_ext[jphc_ext.size() - 1][0],
                                    jphc_ext[jphc_ext.size() - 1][1]));
    }
    if (!ext_ring.is_simple()) {
      std::cerr << "ERROR: exterior ring not a simple polygon" << std::endl;
      _Exit(13);
    }

    // We adopt the convention that exterior rings are counterclockwise
    // oriented, interior rings clockwise oriented. If the orientation in the
    // GeoJSON does not match our convention, we reverse the polygon.
    if (ext_ring.is_clockwise_oriented()) {
      std::cout << "Exterior ring is clockwise" << std::endl;
      ext_ring.reverse_orientation();
    }

    // Store interior ring
    std::vector<Polygon_2<K> > int_ring_v;
    for (int i = 1; i < json_pgn_holes_container.size(); i++) {
      Polygon_2<K> int_ring;
      json jphc_int = json_pgn_holes_container[i];
      for (int j = 0; j < jphc_int.size() - 1; j++) {
        int_ring.push_back(K::Point_2(jphc_int[j][0], jphc_int[j][1]));
      }
      int_ring_v.push_back(int_ring);
      if (jphc_int[0][0] != jphc_int[jphc_int.size() - 1][0] ||
          jphc_int[0][1] != jphc_int[jphc_int.size() - 1][1]) {
        int_ring.push_back(K::Point_2(jphc_int[jphc_int.size() - 1][0],
                                      jphc_int[jphc_int.size() - 1][1]));
      }
      if (!int_ring.is_simple()) {
        std::cerr << "ERROR: interior ring not a simple polygon" << std::endl;
        _Exit(14);
      }
      if (int_ring.is_counterclockwise_oriented()) {
        std::cout << "Interior ring is counterclockwise" << std::endl;
        int_ring.reverse_orientation();
      }
    }
    PolygonWH pgnWH(ext_ring, int_ring_v.begin(), int_ring_v.end());

    gd.push_back(pgnWH);
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
