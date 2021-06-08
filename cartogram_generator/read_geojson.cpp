#include "geo_div.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include <nlohmann/json.hpp>
#include "csv.hpp"
#include <iostream>
#include <fstream>

void check_geojson_validity(const nlohmann::json j)
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
  const nlohmann::json features = j["features"];
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
    const nlohmann::json geometry = feature["geometry"];
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

GeoDiv json_to_cgal(const std::string id,
                    const nlohmann::json json_coords_raw,
                    bool is_polygon)
{
  GeoDiv gd(id);
  nlohmann::json json_coords;
  if (is_polygon) {
    json_coords["0"] = json_coords_raw;
  } else {
    json_coords = json_coords_raw;
  }
  for (auto json_pgn_holes_container : json_coords) {

    // Store exterior ring in CGAL format
    Polygon ext_ring;
    const nlohmann::json jphc_ext = json_pgn_holes_container[0];
    for (unsigned int j = 0; j < jphc_ext.size() - 1; ++j) {
      ext_ring.push_back(Point((double)jphc_ext[j][0],
                               (double)jphc_ext[j][1]));
    }

    // CGAL considers a polygon as simple only if first vertex and last vertex
    // are different
    const unsigned int last_index = jphc_ext.size() - 1;
    if (jphc_ext[0][0] != jphc_ext[last_index][0] ||
        jphc_ext[0][1] != jphc_ext[last_index][1]) {
      ext_ring.push_back(Point((double)jphc_ext[last_index][0],
                               (double)jphc_ext[last_index][1]));
    }
    if (!ext_ring.is_simple()) {
      std::cerr << "ERROR: exterior ring not a simple polygon" << std::endl;
      _Exit(13);
    }

    // We adopt the convention that exterior rings are counterclockwise
    // oriented, interior rings clockwise oriented. If the orientation in the
    // GeoJSON does not match our convention, we reverse the polygon.
    if (ext_ring.is_clockwise_oriented()) {
      ext_ring.reverse_orientation();
    }

    // Store interior ring
    std::vector<Polygon> int_ring_v;
    for (unsigned int i = 1; i < json_pgn_holes_container.size(); ++i) {
      Polygon int_ring;
      const nlohmann::json jphc_int = json_pgn_holes_container[i];
      for (unsigned int j = 0; j < jphc_int.size() - 1; ++j) {
        int_ring.push_back(Point((double)jphc_int[j][0],
                                 (double)jphc_int[j][1]));
      }
      const unsigned int last_index = jphc_int.size() - 1;
      if (jphc_int[0][0] != jphc_int[last_index][0] ||
          jphc_int[0][1] != jphc_int[last_index][1]) {
        int_ring.push_back(Point((double)jphc_int[last_index][0],
                                 (double)jphc_int[last_index][1]));
      }
      if (!int_ring.is_simple()) {
        std::cerr << "ERROR: interior ring not a simple polygon" << std::endl;
        _Exit(14);
      }
      if (int_ring.is_counterclockwise_oriented()) {
        int_ring.reverse_orientation();
      }
      int_ring_v.push_back(int_ring);
    }
    const Polygon_with_holes pwh(ext_ring,
                                 int_ring_v.begin(),
                                 int_ring_v.end());
    gd.push_back(pwh);
  }
  return gd;
}

void print_properties_map(std::map<std::string, std::vector<std::string>>
                          properties_map, unsigned long chosen_number)
{
  unsigned long i = 0;
  for (auto [key, value_vec] : properties_map) {
    i++;
    if (chosen_number == i || chosen_number == properties_map.size() + 1) {
      std::cout << i << ". " << key << ": { ";
      for (long unsigned int j = 0; j < value_vec.size(); j++) {
          std::cout << value_vec[j];
        if (j < value_vec.size() - 1 && j < 5) {
          std::cout << ", ";
        } else {
          if (j < value_vec.size() - 1) {
            std::cout << " ...";
          } else {
            std::cout << " }";
          }
          break;
        }
      }
      std::cout << std::endl;
    }
  }
}

void read_geojson(const std::string geometry_file_name,
                  CartogramInfo *cart_info,
                  bool make_csv)
{
  bool is_polygon;
  bool polygon_warning_has_been_issued = false;

  // Open file
  std::ifstream in_file(geometry_file_name);
  if (!in_file) {
    throw std::system_error(errno,
                            std::system_category(),
                            "failed to open " + geometry_file_name);
  }

  // Parse JSON
  nlohmann::json j;
  try {
    in_file >> j;
  } catch (nlohmann::json::parse_error& e) {
    std::cerr << "ERROR: " << e.what() << '\n'
              << "exception id: " << e.id << '\n'
              << "byte position of error: " << e.byte << std::endl;
    _Exit(3);
  }
  check_geojson_validity(j);
  std::set<std::string> ids_in_geojson;

  // Iterate through each inset
  for (auto &inset_state : *cart_info->ref_to_inset_states()) {
    for (auto feature : j["features"]) {
      const nlohmann::json geometry = feature["geometry"];
      is_polygon = (geometry["type"] == "Polygon");
      if (is_polygon && !polygon_warning_has_been_issued) {
        std::cout << "Warning: support for Polygon geometry experimental, "
                  << "for best results use MultiPolygon" << "\n";
        polygon_warning_has_been_issued = true;
      }
      if (!make_csv) {

        // Storing ID from properties
        const nlohmann::json properties = feature["properties"];
        if (!properties.contains(cart_info->id_header()) &&
            cart_info->id_header() != "") { // Visual file not provided
          std::cerr << "ERROR: In GeoJSON, there is no property "
                    << cart_info->id_header()
                    << " in feature." << std::endl;
          std::cerr << "Available properties are: "
                    << properties
                    << std::endl;
          _Exit(16);
        }

        // Use dump() instead of get() so that we can handle string and numeric
        // IDs in GeoJSON. Both types of IDs are converted to C++ strings.
        std::string id = properties[cart_info->id_header()].dump();

        // We only need to check whether the front of the string is '"' because
        // dump automatically prefixes and postfixes a '"' to any non-NULL string
        // that is not an integer
        if (id.front() == '"') {
          id = id.substr(1, id.length() - 2);
        }
        if (inset_state.pos() == cart_info->inset_at_gd(id)) {
          if (ids_in_geojson.contains(id)) {
            std::cerr << "ERROR: ID "
                      << id
                      << " appears more than once in GeoJSON"
                      << std::endl;
            _Exit(17);
          }
          if (id == "null") {
            std::cerr << "ERROR: ID in GeoJSON is null" << std::endl;
            _Exit(18);
          }
          ids_in_geojson.insert(id);
          const GeoDiv gd = json_to_cgal(id, geometry["coordinates"], is_polygon);
          inset_state.push_back(gd);
        }
      }
    }
  }

  // Creating a CSV from the given GeoJSON file
  if (make_csv) {

    // Declare map for key-value pairs
    std::map<std::string, std::vector<std::string>> properties_map;
    for (auto feature : j["features"]) {
      for (auto property_item : feature["properties"].items()) {
        auto key = property_item.key();

        // Handling strings, and numbers
        std::string value = property_item.value().dump();
        if (value.front() == '"') {
          value = value.substr(1, value.length() - 2);
        }
        auto value_vec = properties_map[key];
        bool value_not_inside =
          std::find(value_vec.begin(), value_vec.end(), value)
          == value_vec.end();
        if (value != "null" && !value.empty() && value_not_inside) {
          properties_map[key].push_back(value);
        }
      }
    }

    // Discard keys with repeating or missing values
    std::map<std::string, std::vector<std::string>>
        viable_properties_map = properties_map;
    for (auto [key, value_vec] : properties_map) {
      if (value_vec.size() < j["features"].size()) {
        viable_properties_map.erase(key);
      }
    }
    std::cout << std::endl;

    // Present user with all possible identifiers and a few examples
    std::cout << "These are the unique identifiers and their values: ";
    std::cout << std::endl;
    print_properties_map(
      viable_properties_map, viable_properties_map.size() + 1);
    std::cout << viable_properties_map.size() + 1 << ". All";
    std::cout << std::endl << std::endl;

    // Have the user choose which key(s) they want to use as the identifier(s)
    std::cout << "Please enter your number here: ";
    unsigned long chosen_number = 0;
    while (std::cin.fail()
           || chosen_number < 1
           || chosen_number > viable_properties_map.size() + 1) {

      // Prompting User for Input
      std::cout << "Please enter your number here: ";
      std::cin >> chosen_number;
      if (std::cin.fail()) {
        std::cout << "Invalid input! Try again." << std::endl;

        // Clearing std::cin buffer
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      } else if (chosen_number < 1
                 || chosen_number > viable_properties_map.size() + 1) {
        std::cout << "Please enter a number between 1 and "
                  << viable_properties_map.size() + 1
                  << std::endl;
      }
    }
    std::cout << std::endl;

    // Declare chosen identifier(s)
    std::map<std::string, std::vector<std::string>> chosen_identifiers;
    size_t i = 0;
    for (auto [key, value_vec] : viable_properties_map) {
      i++;
      if (chosen_number == i
          || chosen_number == viable_properties_map.size() + 1) {
        chosen_identifiers[key] = value_vec;
      }
    }

    // Printing chosen identifiers
    std::cout << "Chosen identifiers: " << std::endl;
    print_properties_map(viable_properties_map, chosen_number);
    std::cout << std::endl;

    // Writing CSV
    std::ofstream out_file_csv;
    out_file_csv.open ("template_from_geojson.csv");
    if (!out_file_csv) {
      throw std::system_error(errno,
                              std::system_category(),
                              "failed to open template_from_geojson.csv");
    }

    // Each vector of strings will represent one row
    std::vector<std::vector<std::string>>
      csv_rows(chosen_identifiers.begin()->second.size() + 1);

    // Converting map into a vector
    int column = 0;
    for (auto [column_name, ids] : chosen_identifiers) {
      csv_rows[0].push_back(column_name);
      if (column == 0) {
        csv_rows[0].push_back("Cartogram Data (eg. Population)");
        csv_rows[0].push_back("Color");
        csv_rows[0].push_back("Inset");
        csv_rows[0].push_back("Abbreviation");
      }
      for (size_t k = 0; k < ids.size(); k++) {
        csv_rows[k + 1].push_back(ids[k]);
        if (column == 0) {
          for (size_t i = 0; i < 4; i++) {
            csv_rows[k + 1].push_back("");
          }
        }
      }
      column++;
    }

    // Writing to CSV writer object
    auto writer = csv::make_csv_writer(out_file_csv);
    for (auto row : csv_rows) {
      writer << row;
    }

    // Closing out_file and exiting
    out_file_csv.close();
    _Exit(19);
  }

  // Check whether all IDs in visual_variable_file appear in GeoJSON
  const std::set<std::string> ids_in_vv_file =
    cart_info->ids_in_visual_variables_file();
  std::set<std::string> ids_not_in_geojson;
  std::set_difference(ids_in_vv_file.begin(), ids_in_vv_file.end(),
                      ids_in_geojson.begin(), ids_in_geojson.end(),
                      std::inserter(ids_not_in_geojson,
                                    ids_not_in_geojson.end()));
  if (!ids_not_in_geojson.empty()) {
    std::cerr << "ERROR: Mismatch between GeoJSON and "
              << cart_info->visual_variable_file()
              << "."
              << std::endl;
    std::cerr << "The following IDs do not appear in the GeoJSON or CSV:"
              << std::endl;
    for (auto id : ids_not_in_geojson) {
      std::cerr << "  " << id << std::endl;
    }
    _Exit(20);
  }

  // Check whether all IDs in GeoJSON appear in visual_variable_file.
  std::set<std::string> ids_not_in_vv;
  std::set_difference(ids_in_geojson.begin(), ids_in_geojson.end(),
                      ids_in_vv_file.begin(), ids_in_vv_file.end(),
                      std::inserter(ids_not_in_vv,
                                    ids_not_in_vv.end()));
  if (!ids_not_in_vv.empty()) {
    std::cerr << "ERROR: Mismatch between GeoJSON and "
              << cart_info->visual_variable_file()
              << "."
              << std::endl;
    std::cerr << "The following IDs do not appear in "
              << cart_info->visual_variable_file()
              << ": "
              << std::endl;
    for (auto id : ids_not_in_vv) {
      std::cerr << "  " << id << std::endl;
    }
    _Exit(21);
  }
  return;
}
