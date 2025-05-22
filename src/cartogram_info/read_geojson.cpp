#include "cartogram_info.hpp"
#include "constants.hpp"
#include "csv.hpp"

inline std::string strip_quotes(const std::string &s)
{
  if (s.front() == '"' && s.back() == '"') {
    return s.substr(1, s.size() - 2);
  }
  return s;
}

static void check_geojson_validity(const nlohmann::json &j)
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
  const auto &features = j["features"];
  for (const auto &feature : features) {
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
    const auto &geometry = feature["geometry"];
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
}

static std::pair<GeoDiv, bool> json_to_geodiv(
  const std::string &id,
  const nlohmann::json &json_coords_raw,
  const bool is_polygon)
{
  GeoDiv gd(id);
  nlohmann::json json_coords;
  if (is_polygon) {
    json_coords["0"] = json_coords_raw;
  } else {
    json_coords = json_coords_raw;
  }
  bool erico;  // Exterior ring is clockwise oriented?
  for (const auto &json_pgn_holes_container : json_coords) {

    // Store exterior ring in CGAL format
    Polygon ext_ring;
    const auto jphc_ext = json_pgn_holes_container[0];
    for (size_t j = 0; j < jphc_ext.size() - 1; ++j) {
      ext_ring.push_back(Point(
        static_cast<double>(jphc_ext[j][0]),
        static_cast<double>(jphc_ext[j][1])));
    }

    // CGAL considers a polygon as simple only if first vertex and last vertex
    // are different
    const auto last_ext_index = jphc_ext.size() - 1;
    if (
      jphc_ext[0][0] != jphc_ext[last_ext_index][0] ||
      jphc_ext[0][1] != jphc_ext[last_ext_index][1]) {
      ext_ring.push_back(Point(
        static_cast<double>(jphc_ext[last_ext_index][0]),
        static_cast<double>(jphc_ext[last_ext_index][1])));
    }
    if (!ext_ring.is_simple()) {
      std::cerr << "ERROR: (GeoJSON Parsing) exterior ring not a simple polygon" << std::endl;
      _Exit(13);
    }

    // We adopt the convention that exterior rings are counterclockwise
    // oriented, interior rings clockwise oriented. If the orientation in the
    // GeoJSON does not match our convention, we reverse the polygon.

    // TODO: Currently, only the last exterior ring in the GeoJSON determines
    //       whether original_ext_ring_is_clockwise in CartogramInfo is true.
    //       This strategy works for most geospatial boundary files in the
    //       wild, but it would still be sensible to allow cases where there
    //       are external rings with opposite winding directions.
    erico = ext_ring.is_clockwise_oriented();
    if (erico) {
      ext_ring.reverse_orientation();
    }

    // Store interior ring
    std::vector<Polygon> int_ring_v;
    for (size_t i = 1; i < json_pgn_holes_container.size(); ++i) {
      Polygon int_ring;
      const auto jphc_int = json_pgn_holes_container[i];
      for (size_t j = 0; j < jphc_int.size() - 1; ++j) {
        int_ring.push_back(Point(
          static_cast<double>(jphc_int[j][0]),
          static_cast<double>(jphc_int[j][1])));
      }
      const size_t last_int_index = jphc_int.size() - 1;
      if (
        jphc_int[0][0] != jphc_int[last_int_index][0] ||
        jphc_int[0][1] != jphc_int[last_int_index][1]) {
        int_ring.push_back(Point(
          static_cast<double>(jphc_int[last_int_index][0]),
          static_cast<double>(jphc_int[last_int_index][1])));
      }
      if (!int_ring.is_simple()) {
        std::cerr << "ERROR: (GeoJSON Parsing) interior ring not a simple polygon" << std::endl;
        _Exit(14);
      }
      if (int_ring.is_counterclockwise_oriented()) {
        int_ring.reverse_orientation();
      }
      int_ring_v.push_back(int_ring);
    }
    const Polygon_with_holes pwh(
      ext_ring,
      int_ring_v.begin(),
      int_ring_v.end());
    gd.push_back(pwh);
  }
  return {gd, erico};
}

static void print_properties_map(
  const std::map<std::string, std::vector<std::string>> &properties_map,
  const unsigned long chosen_number)
{
  const size_t max_n_printed_values = 5;
  const auto value_vec = properties_map.begin()->second;
  const size_t n_printed_values = std::min(
    value_vec.size(),
    static_cast<unsigned long>(max_n_printed_values));
  size_t i = 0;
  for (const auto &[key, val] : properties_map) {
    ++i;
    if (chosen_number == i || chosen_number == properties_map.size() + 1) {
      std::cerr << i << ". " << key << ": {";
      for (unsigned long j = 0; j < n_printed_values - 1; ++j) {
        std::cerr << val[j] << ", ";
      }
      std::cerr << val[n_printed_values - 1];
      if (val.size() > n_printed_values) {
        std::cerr << " ...";
      }
      std::cerr << "}" << std::endl;
    }
  }
}

static nlohmann::json load_geojson(const std::string &geometry_file_name)
{
  std::ifstream in_file(geometry_file_name);
  if (!in_file) {
    std::cerr << "ERROR reading GeoJSON: failed to open " << geometry_file_name
              << std::endl;
    _Exit(3);
  }

  nlohmann::json j;
  try {
    in_file >> j;
  } catch (nlohmann::json::parse_error &e) {
    std::cerr << "ERROR: " << e.what() << ". Exception id: " << e.id
              << ". Byte position of error: " << e.byte << std::endl;
    _Exit(3);
  }

  return j;
}

// Read coordinate reference system if it is included in the GeoJSON
static void extract_crs(const nlohmann::json &j, std::string &crs)
{
  if (
    j.contains("crs") && j["crs"].contains("properties") &&
    j["crs"]["properties"].contains("name")) {
    crs = j["crs"]["properties"]["name"];
  }
  std::cerr << "Coordinate reference system: " << crs << std::endl;
}

static std::map<std::string, std::vector<std::string>> extract_unique_properties_map(
  const nlohmann::json &j)
{
  std::map<std::string, std::vector<std::string>> properties_map;
  for (const auto &feature : j["features"]) {
    for (const auto &property_item : feature["properties"].items()) {
      const std::string key = property_item.key();
      const std::string value = strip_quotes(property_item.value().dump());

      const std::vector<std::string> value_vec = properties_map[key];
      const bool value_not_inside =
        std::find(value_vec.begin(), value_vec.end(), value) ==
        value_vec.end();
      if (value != "null" && !value.empty() && value_not_inside)
        properties_map[key].push_back(value);
    }
  }

  // Discard keys with repeating or missing values
  auto unique_properties_map = properties_map;
  for (const auto &[key, value_vec] : properties_map) {
    if (value_vec.size() < j["features"].size())
      unique_properties_map.erase(key);
  }

  if (unique_properties_map.empty()) {
    std::cerr << "ERROR: No unique properties found in GeoJSON" << std::endl;
    _Exit(15);
  }

  return unique_properties_map;
}

static void generate_csv_template(
  const nlohmann::json &j,
  const std::string &map_name_)
{
  std::map<std::string, std::vector<std::string>> viable_properties_map =
    extract_unique_properties_map(j);

  std::cerr << std::endl;

  // Have the users choose which key(s) they want to use as the
  // identifier(s) if more than one key available
  unsigned long chosen_number = 0;
  if (viable_properties_map.size() > 1) {

    // Present user with all possible identifiers and a few examples
    std::cerr << "These are the unique identifiers and their values:\n"
              << std::endl;
    print_properties_map(
      viable_properties_map,
      viable_properties_map.size() + 1);
    std::cerr << viable_properties_map.size() + 1 << ". All\n" << std::endl;
    while (std::cin.fail() || chosen_number < 1 ||
           chosen_number > viable_properties_map.size() + 1) {

      // Prompt user for input
      std::cerr << "Please enter your number here: ";
      std::cin >> chosen_number;
      if (std::cin.fail()) {
        std::cerr << "Invalid input! Try again." << std::endl;

        // Clear std::cin buffer
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      } else if (
        chosen_number < 1 ||
        chosen_number > viable_properties_map.size() + 1) {
        std::cerr << "Please enter a number between 1 and "
                  << viable_properties_map.size() + 1 << std::endl;
      }
    }
    std::cerr << std::endl;
  } else {
    std::cerr << "Only one unique identifier found: ";
    print_properties_map(
      viable_properties_map,
      viable_properties_map.size() + 1);
    std::cerr << std::endl;
    ++chosen_number;
  }

  // Declare chosen identifier(s)
  std::map<std::string, std::vector<std::string>> chosen_identifiers;
  size_t i = 0;
  for (const auto &[key, value_vec] : viable_properties_map) {
    ++i;
    if (
      chosen_number == i ||
      chosen_number == viable_properties_map.size() + 1) {
      chosen_identifiers[key] = value_vec;
    }
  }

  // Print chosen identifier(s)
  std::cerr << "Chosen identifier(s): " << std::endl;
  print_properties_map(viable_properties_map, chosen_number);
  std::cerr << std::endl;

  // Write CSV
  std::ofstream out_file_csv;
  const auto csv_name = map_name_ + ".csv";
  out_file_csv.open(csv_name);
  if (!out_file_csv) {
    std::cerr
      << "ERROR writing GeoJSON: failed to open template_from_geojson.csv"
      << std::endl;
  }

  // Each vector of strings will represent one row
  std::vector<std::vector<std::string>> csv_rows(
    chosen_identifiers.begin()->second.size() + 1);

  // Converting map into a vector
  size_t column = 0;
  for (const auto &[column_name, ids] : chosen_identifiers) {
    csv_rows[0].push_back(column_name);
    if (column == 0) {
      csv_rows[0].push_back("Cartogram Data (eg. Population)");
      csv_rows[0].push_back("Color");
      csv_rows[0].push_back("Inset");
      csv_rows[0].push_back("Label");
    }
    for (size_t k = 0; k < ids.size(); ++k) {
      csv_rows[k + 1].push_back(ids[k]);
      if (column == 0) {
        for (size_t l = 0; l < 4; ++l) {
          csv_rows[k + 1].push_back("");
        }
      }
    }
    ++column;
  }

  // Write to CSV object
  auto writer = csv::make_csv_writer(out_file_csv);
  for (const auto &row : csv_rows) {
    writer << row;
  }

  // Close out_file and exit
  out_file_csv.close();
}

static std::vector<std::string> extract_initial_order_of_ids(
  const nlohmann::json &j,
  const std::string &id_header)
{
  std::vector<std::string> initial_id_order;
  for (const auto &feature : j["features"]) {
    const auto properties = feature["properties"];
    assert(properties.contains(id_header));
    const auto id = strip_quotes(properties[id_header].dump());
    initial_id_order.push_back(id);
  }
  return initial_id_order;
}

void CartogramInfo::construct_inset_state_from_geodivs(const nlohmann::json &j)
{
  InsetState inset_state("C", args_);
  for (const auto &feature : j["features"]) {
    const auto geometry = feature["geometry"];
    const bool is_polygon = (geometry["type"] == "Polygon");

    std::string id = strip_quotes(feature["properties"][id_header_].dump());
    auto [gd, erico] = json_to_geodiv(id, geometry["coordinates"], is_polygon);
    inset_state.push_back(gd);
    gd_to_inset_.emplace(id, "C");
    original_ext_ring_is_clockwise_ = erico;
  }
  inset_states_.emplace_back(std::move(inset_state));
}

void CartogramInfo::read_geojson()
{
  std::string geometry_file_name = args_.geo_file_name;
  nlohmann::json j = load_geojson(geometry_file_name);
  check_geojson_validity(j);

  if (args_.make_csv) {
    // Update map_name to be based on geometry file name instead
    set_map_name(geometry_file_name);
    generate_csv_template(j, map_name_);
    _Exit(19);
  }

  extract_crs(j, crs_);
  // Skip projection, this is an output from our program
  if (crs_ == custom_crs) {
    std::cerr << "WARNING: " << custom_crs << " detected. "
              << "Applying --skip_projection flag." << std::endl;
    args_.skip_projection = true;
  }

  unique_properties_map_ = extract_unique_properties_map(j);
  assert(unique_properties_map_.size() > 0);

  // Set the first key inside the unique_properties_map_ as the default ID header
  id_header_ = unique_properties_map_.begin()->first;

  initial_id_order_ = extract_initial_order_of_ids(j, id_header_);

  construct_inset_state_from_geodivs(j);
}
