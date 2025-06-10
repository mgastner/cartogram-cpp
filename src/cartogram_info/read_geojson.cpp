#include "cartogram_info.hpp"
#include "constants.hpp"
#include "csv.hpp"
#include "round_point.hpp"
#include "simdjson.h"

inline bool has_key(
  const simdjson::dom::object &obj,
  std::string_view key) noexcept
{
  return obj.at_key(key).error() == simdjson::SUCCESS;
}

inline simdjson::dom::element at_key(
  const simdjson::dom::object &obj,
  std::string_view key)
{
  auto e = obj.at_key(key);
  if (e.error() != simdjson::SUCCESS) {
    throw simdjson::simdjson_error(e.error());
  }
  return e;
}

inline std::string to_std_string(const simdjson::dom::element &e)
{
  if (e.is_string()) {
    auto sv = e.get_string();
    if (sv.error())
      throw simdjson::simdjson_error(sv.error());
    return std::string{sv.value()};
  }
  if (e.is_null())
    return "null";
  return simdjson::to_string(e);
}

inline double as_double(const simdjson::dom::element &e)
{
  return e.get_double().value();
}

inline simdjson::dom::array as_array(const simdjson::dom::element &e)
{
  return e.get_array().value();
}

inline std::string strip_quotes(const std::string &s)
{
  if (s.front() == '"' && s.back() == '"') {
    return s.substr(1, s.size() - 2);
  }
  return s;
}

static void check_geojson_validity(const simdjson::dom::element &j)
{
  const simdjson::dom::object obj = j.get_object();

  if (!has_key(obj, "type")) {
    std::cerr << "ERROR: JSON does not contain a key 'type'" << std::endl;
    std::exit(4);
  }

  if (to_std_string(at_key(obj, "type")) != "FeatureCollection") {
    std::cerr << "ERROR: JSON is not a valid GeoJSON FeatureCollection"
              << std::endl;
    std::exit(5);
  }

  if (!has_key(obj, "features")) {
    std::cerr << "ERROR: JSON does not contain a key 'features'" << std::endl;
    std::exit(6);
  }

  const simdjson::dom::array features = at_key(obj, "features");
  for (const simdjson::dom::element &feature_elem : features) {
    const simdjson::dom::object feature = feature_elem.get_object();

    if (!has_key(feature, "type")) {
      std::cerr << "ERROR: JSON contains a 'Features' element without key "
                << "'type'" << std::endl;
      std::exit(7);
    }

    if (to_std_string(at_key(feature, "type")) != "Feature") {
      std::cerr << "ERROR: JSON contains a 'Features' element whose type "
                << "is not 'Feature'" << std::endl;
      std::exit(8);
    }

    if (!has_key(feature, "geometry")) {
      std::cerr << "ERROR: JSON contains a feature without key 'geometry'"
                << std::endl;
      std::exit(9);
    }

    const simdjson::dom::object geometry = at_key(feature, "geometry");
    if (!has_key(geometry, "type")) {
      std::cerr << "ERROR: JSON contains geometry without key 'type'"
                << std::endl;
      std::exit(10);
    }

    if (!has_key(geometry, "coordinates")) {
      std::cerr << "ERROR: JSON contains geometry without key 'coordinates'"
                << std::endl;
      std::exit(11);
    }

    const std::string gtype = to_std_string(at_key(geometry, "type"));
    if (gtype != "MultiPolygon" && gtype != "Polygon") {
      std::cerr << "ERROR: JSON contains unsupported geometry " << gtype
                << std::endl;
      std::exit(12);
    }
  }
}

static std::pair<GeoDiv, bool> json_to_geodiv(
  const std::string &id,
  const simdjson::dom::element &json_coords_raw,
  const bool is_polygon)
{
  GeoDiv gd(id);
  bool erico = false;  // Exterior ring is clockwise oriented?

  auto process_polygon = [&](const simdjson::dom::array &rings) {
    // Store exterior ring in CGAL format
    const simdjson::dom::array ext = as_array(rings.at(0));
    Polygon ext_ring;

    for (size_t j = 0; j + 1 < ext.size(); ++j) {
      const simdjson::dom::array pt = as_array(ext.at(j));
      ext_ring.push_back(Point(as_double(pt.at(0)), as_double(pt.at(1))));
    }

    // CGAL considers a polygon as simple only if first vertex and last vertex
    // are different
    const size_t last_ext_index = ext.size() - 1;
    const simdjson::dom::array pt0 = as_array(ext.at(0));
    const simdjson::dom::array lastp = as_array(ext.at(last_ext_index));

    if (
      !almost_equal(as_double(pt0.at(0)), as_double(lastp.at(0))) ||
      !almost_equal(as_double(pt0.at(1)), as_double(lastp.at(1)))) {
      ext_ring.push_back(
        Point(as_double(lastp.at(0)), as_double(lastp.at(1))));
    }

    if (!ext_ring.is_simple()) {
      std::cerr
        << "ERROR: (GeoJSON Parsing) exterior ring not a simple polygon"
        << std::endl;
      std::exit(13);
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
    if (erico)
      ext_ring.reverse_orientation();

    // Store interior ring
    std::vector<Polygon> int_ring_v;
    for (size_t i = 1; i < rings.size(); ++i) {
      const simdjson::dom::array int_r = as_array(rings.at(i));
      Polygon int_ring;

      for (size_t j = 0; j + 1 < int_r.size(); ++j) {
        const simdjson::dom::array pt = as_array(int_r.at(j));
        int_ring.push_back(Point(as_double(pt.at(0)), as_double(pt.at(1))));
      }

      const size_t last_int_index = int_r.size() - 1;
      const simdjson::dom::array pt0_i = as_array(int_r.at(0));
      const simdjson::dom::array lastp_i = as_array(int_r.at(last_int_index));

      if (
        !almost_equal(as_double(pt0_i.at(0)), as_double(lastp_i.at(0))) ||
        !almost_equal(as_double(pt0_i.at(1)), as_double(lastp_i.at(1)))) {
        int_ring.push_back(
          Point(as_double(lastp_i.at(0)), as_double(lastp_i.at(1))));
      }

      if (!int_ring.is_simple()) {
        std::cerr
          << "ERROR: (GeoJSON Parsing) interior ring not a simple polygon"
          << std::endl;
        std::exit(14);
      }

      if (int_ring.is_counterclockwise_oriented())
        int_ring.reverse_orientation();
      int_ring_v.push_back(int_ring);
    }

    gd.push_back(
      Polygon_with_holes(ext_ring, int_ring_v.begin(), int_ring_v.end()));
  };

  if (is_polygon) {
    process_polygon(as_array(json_coords_raw));
  } else {
    const simdjson::dom::array polys = as_array(json_coords_raw);
    for (const simdjson::dom::element &poly_elem : polys)
      process_polygon(as_array(poly_elem));
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

static simdjson::dom::element load_geojson(
  const std::string &geometry_file_name)
{
  static simdjson::dom::parser parser;
  static std::vector<simdjson::padded_string> buffers;

  try {
    auto &buf =
      buffers.emplace_back(simdjson::padded_string::load(geometry_file_name));
    return parser.parse(buf);
  } catch (const simdjson::simdjson_error &e) {
    std::cerr << "ERROR reading GeoJSON: " << e.what() << std::endl;
    std::exit(3);
  }
}

// Read coordinate reference system if it is included in the GeoJSON
static void extract_crs(const simdjson::dom::element &j, std::string &crs)
{
  if (j.type() != simdjson::dom::element_type::OBJECT)
    return;

  const simdjson::dom::object obj = j.get_object();
  if (has_key(obj, "crs")) {
    const simdjson::dom::object crs_obj = at_key(obj, "crs");
    if (has_key(crs_obj, "properties")) {
      const simdjson::dom::object props = at_key(crs_obj, "properties");
      if (has_key(props, "name")) {
        crs = to_std_string(at_key(props, "name"));
      }
    }
  }
  std::cerr << "Coordinate reference system: " << crs << std::endl;
}

static std::map<std::string, std::vector<std::string>>
extract_unique_properties_map(const simdjson::dom::element &j)
{
  std::map<std::string, std::vector<std::string>> properties_map;

  const simdjson::dom::array features = at_key(j.get_object(), "features");
  for (const simdjson::dom::element &feature_elem : features) {
    const simdjson::dom::object props =
      at_key(feature_elem.get_object(), "properties");

    for (auto [k, v] : props) {
      const std::string key = std::string{k};
      const std::string value = to_std_string(v);

      const std::vector<std::string> &value_vec = properties_map[key];
      const bool value_not_inside =
        std::find(value_vec.begin(), value_vec.end(), value) ==
        value_vec.end();
      if (value != "null" && !value.empty() && value_not_inside)
        properties_map[key].push_back(value);
    }
  }

  // Discard keys with repeating or missing values
  auto unique_properties_map = properties_map;
  const size_t n_features = features.size();
  for (const auto &[key, value_vec] : properties_map) {
    if (value_vec.size() < n_features)
      unique_properties_map.erase(key);
  }

  if (unique_properties_map.empty()) {
    std::cerr << "ERROR: No unique properties found in GeoJSON" << std::endl;
    std::exit(15);
  }

  return unique_properties_map;
}

static void generate_csv_template(
  const simdjson::dom::element &j,
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
  const simdjson::dom::element &j,
  const std::string &id_header)
{
  std::vector<std::string> initial_id_order;

  const simdjson::dom::array features = at_key(j.get_object(), "features");
  for (const simdjson::dom::element &feature_elem : features) {
    const simdjson::dom::object properties =
      at_key(feature_elem.get_object(), "properties");
    assert(has_key(properties, id_header));
    const auto id = strip_quotes(to_std_string(at_key(properties, id_header)));
    initial_id_order.push_back(id);
  }
  return initial_id_order;
}

void CartogramInfo::construct_inset_state_from_geodivs(
  const simdjson::dom::element &j)
{
  InsetState inset_state("C", args_);
  const simdjson::dom::array features = at_key(j.get_object(), "features");

  for (const simdjson::dom::element &feature_elem : features) {
    const simdjson::dom::object geometry =
      at_key(feature_elem.get_object(), "geometry");
    const bool is_polygon =
      (to_std_string(at_key(geometry, "type")) == "Polygon");

    std::string id = strip_quotes(to_std_string(
      at_key(at_key(feature_elem.get_object(), "properties"), id_header_)));

    auto [gd, erico] =
      json_to_geodiv(id, at_key(geometry, "coordinates"), is_polygon);
    inset_state.push_back(gd);
    gd_to_inset_.emplace(id, "C");
    original_ext_ring_is_clockwise_ = erico;
  }
  inset_states_.emplace_back(std::move(inset_state));
}

void CartogramInfo::read_geojson()
{
  std::string geometry_file_name = args_.geo_file_name;
  simdjson::dom::element j = load_geojson(geometry_file_name);
  check_geojson_validity(j);

  if (args_.make_csv) {
    // Update map_name to be based on geometry file name instead
    set_map_name(geometry_file_name);
    generate_csv_template(j, map_name_);
    std::exit(19);
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

  // Set the first key inside the unique_properties_map_ as the default ID
  // header
  id_header_ = unique_properties_map_.begin()->first;

  initial_id_order_ = extract_initial_order_of_ids(j, id_header_);

  construct_inset_state_from_geodivs(j);
}
