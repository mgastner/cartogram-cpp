#include "cartogram_info.hpp"
#include "string_to_decimal_converter.hpp"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <jsoncons/json.hpp>
#include <jsoncons_ext/csv/csv.hpp>
#include <locale>
#include <map>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

static inline void ltrim_ascii(std::string &s)
{
  size_t i = 0;
  while (i < s.size()) {
    unsigned char ch = static_cast<unsigned char>(s[i]);
    if (
      ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || ch == '\v' ||
      ch == '\f')
      ++i;
    else
      break;
  }
  if (i)
    s.erase(0, i);
}

static inline void rtrim_ascii(std::string &s)
{
  while (!s.empty()) {
    unsigned char ch = static_cast<unsigned char>(s.back());
    if (
      ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || ch == '\v' ||
      ch == '\f')
      s.pop_back();
    else
      break;
  }
}

static inline std::string normalize_id_token(std::string s)
{
  // Strip UTF-8 BOM if present
  if (
    s.size() >= 3 && static_cast<unsigned char>(s[0]) == 0xEF &&
    static_cast<unsigned char>(s[1]) == 0xBB &&
    static_cast<unsigned char>(s[2]) == 0xBF) {
    s.erase(0, 3);
  }
  ltrim_ascii(s);
  rtrim_ascii(s);
  return s;
}

static void check_validity_of_area_str(const std::string &area_as_str)
{
  std::string area_process_str = area_as_str.empty() ? "NA" : area_as_str;

  if (!StringToDecimalConverter::is_str_valid_characters(area_process_str)) {
    std::cerr
      << "ERROR: Invalid area string: " << area_process_str
      << ". Area string must only contain 0-9, '.', '-' and ',' or 'NA'."
      << std::endl;
    std::exit(18);
  }

  if (
    !StringToDecimalConverter::is_str_NA(area_process_str) &&
    !StringToDecimalConverter::is_str_correct_format(area_process_str)) {
    std::cerr << "ERROR: Invalid area string format: " << area_process_str
              << std::endl;
    std::exit(19);
  }

  if (area_process_str.front() == '-') {
    std::cerr << "ERROR: Negative area in CSV" << std::endl;
    std::exit(101);
  }
}

static std::string process_inset_pos_str(const std::string &inset_pos_as_str)
{
  std::string inset_pos = inset_pos_as_str.empty() ? "C" : inset_pos_as_str;
  inset_pos[0] =
    static_cast<char>(std::toupper(static_cast<unsigned char>(inset_pos[0])));
  if (inset_pos == "U")
    inset_pos = "T";
  if (inset_pos == "D")
    inset_pos = "B";
  return inset_pos;
}

static void check_validity_of_inset_pos(
  const std::string &inset_pos,
  const std::string &id)
{
  std::unordered_set<std::string> permitted_pos{"C", "L", "R", "T", "B"};
  if (!permitted_pos.contains(inset_pos)) {
    std::cerr << "Unrecognized inset position : " << inset_pos
              << " for Region: " << id << "\nSetting " << id
              << "\'s inset position to Center (C)." << std::endl;
    std::exit(20);
  }
}

using jsoncons::ojson;
namespace jcc = jsoncons::csv;

struct HeaderMap {
  std::vector<std::string> raw;  // as in file
  std::vector<std::string> norm;  // trimmed/BOM-stripped view

  int index_of_norm(const std::string &name) const
  {
    std::string target = name;
    ltrim_ascii(target);
    rtrim_ascii(target);
    // Do not strip BOM from 'name'; caller passes canonical tokens like
    // "Color"
    auto it = std::find(norm.begin(), norm.end(), target);
    return it == norm.end() ? -1 : static_cast<int>(it - norm.begin());
  }

  const std::string &raw_at(size_t i) const
  {
    return raw[i];
  }

  const std::string &norm_at(size_t i) const
  {
    return norm[i];
  }
};

static ojson load_csv_rows(const std::string &path)
{
  std::ifstream is(path);
  if (!is) {
    std::cerr << "ERROR: Cannot open CSV file: " << path << std::endl;
    std::exit(17);
  }
  jcc::csv_options opts;
  opts.assume_header(true).ignore_empty_lines(true).trim(true).infer_types(
    false);  // keep all fields as strings for stable behavior
  try {
    return jcc::decode_csv<ojson>(is, opts);  // array of row-objects
  } catch (const std::exception &e) {
    std::cerr << "ERROR: Failed to parse CSV: " << e.what() << std::endl;
    std::exit(17);
  }
}

static HeaderMap headers_from(const ojson &rows)
{
  HeaderMap hm;
  if (!rows.is_array() || rows.empty())
    return hm;

  const ojson &first = rows.at(0);
  for (const auto &kv : first.object_range()) {
    hm.raw.emplace_back(std::string(kv.key()));
    hm.norm.emplace_back(normalize_id_token(std::string(kv.key())));
  }
  return hm;
}

static int extract_color_col_index(
  const HeaderMap &hm,
  const std::string &color_col_name)
{
  int idx = hm.index_of_norm(color_col_name);
  if (idx < 0 && color_col_name == "Color")
    idx = hm.index_of_norm("Colour");
  return idx;
}

// Find the matching ID columns in both the CSV and GeoJSON file
// Returns the header name of the matching ID column in the GeoJSON file
std::string CartogramInfo::match_id_columns(
  const std::optional<std::string> &id_col)
{
  ojson rows = load_csv_rows(args_.visual_file_name);
  const HeaderMap hm = headers_from(rows);
  std::string csv_id_header_norm;

  // Build normalized sets from GeoJSON unique properties
  std::map<std::string, std::set<std::string>> geojson_properties_info;
  for (auto &[key, properties_vec] : unique_properties_map_) {
    std::set<std::string> s;
    for (const auto &v : properties_vec)
      s.insert(normalize_id_token(v));
    geojson_properties_info[key] = std::move(s);
  }

  std::string matching_id_header;

  auto try_match = [&](const std::string &header_norm) -> bool {
    int idx = hm.index_of_norm(header_norm);
    if (idx < 0)
      return false;
    const std::string &raw_key = hm.raw_at(static_cast<size_t>(idx));

    std::set<std::string> data_set;
    size_t row_count = rows.size();
    for (const auto &r : rows.array_range()) {
      std::string cell;
      if (r.contains(raw_key)) {
        const ojson &v = r.at(raw_key);
        // Values are strings because infer_types(false)
        cell = v.is_string() ? v.as_string() : std::string();
      }
      data_set.insert(normalize_id_token(cell));
    }

    if (data_set.size() != row_count)
      return false;

    for (auto &[key, value_set] : geojson_properties_info) {
      if (data_set == value_set) {
        matching_id_header = key;
        csv_id_header_norm = header_norm;
        return true;
      }
    }
    return false;
  };

  if (id_col) {
    if (!try_match(*id_col)) {
      std::cerr << "Given ID header " << *id_col
                << " does not match with any GeoJSON properties. "
                   "Finding next best matching ID column..."
                << std::endl;
    }
  }

  if (matching_id_header.empty()) {
    for (const auto &h : hm.norm) {
      if (try_match(h)) {
        std::cerr << "Matched ID column: " << h << std::endl;
        break;
      }
    }
  }

  if (matching_id_header.empty()) {
    std::cerr << "ERROR: No valid matching ID header between GeoJSON "
                 "properties and CSV columns could be found."
              << std::endl;
    std::exit(16);
  }

  int idx = hm.index_of_norm(csv_id_header_norm);
  if (idx < 0) {
    std::cerr << "ERROR: Internal error determining ID column index."
              << std::endl;
    std::exit(16);
  }
  id_col_ = idx;  // index into normalized header list
  return matching_id_header;
}

// Updates ID header and inset info
void CartogramInfo::update_id_header_info(
  const std::string &matching_id_header)
{
  std::vector<std::string> old_unique_properties =
    unique_properties_map_[id_header_];
  std::vector<std::string> new_unique_properties =
    unique_properties_map_[matching_id_header];

  std::map<std::string, std::string> geojson_id_to_csv_id;
  for (size_t i = 0; i < old_unique_properties.size(); i++)
    geojson_id_to_csv_id[old_unique_properties[i]] = new_unique_properties[i];

  for (auto &id : initial_id_order_) {
    id = geojson_id_to_csv_id.at(id);
  }

  std::map<std::string, std::string> new_gd_to_inset;
  for (auto &[geojson_id, inset_pos] : gd_to_inset_) {
    const std::string csv_id = geojson_id_to_csv_id.at(geojson_id);
    new_gd_to_inset[csv_id] = inset_pos;
  }
  
  gd_to_inset_ = std::move(new_gd_to_inset);

  for (InsetState &inset_state : inset_states_) {
    inset_state.update_gd_ids(geojson_id_to_csv_id);
  }

  id_header_ = matching_id_header;
}

static void check_validity_of_csv_ids(
  std::map<std::string, std::map<std::string, std::string>> &csv_data,
  const std::vector<std::string> &initial_id_order)
{
  std::vector<std::string> csv_ids;
  csv_ids.reserve(csv_data.size());
  for (const auto &[id, _] : csv_data)
    csv_ids.push_back(id);

  for (const auto &id : csv_ids) {
    if (
      std::find(initial_id_order.begin(), initial_id_order.end(), id) ==
      initial_id_order.end()) {
      std::cerr << "ERROR: ID " << id << " in CSV is not in GeoJSON"
                << std::endl;
      std::exit(21);
    }
  }

  for (const auto &id : initial_id_order) {
    if (std::find(csv_ids.begin(), csv_ids.end(), id) == csv_ids.end()) {
      std::cerr << "WARNING: ID " << id << " in GeoJSON is not in CSV"
                << std::endl;
      csv_data[id] =
        {{"area", "NA"}, {"color", ""}, {"label", ""}, {"inset_pos", "C"}};
    }
  }
}

void CartogramInfo::relocate_geodivs_based_on_inset_pos(
  const std::map<std::string, std::map<std::string, std::string>> &csv_data)
{
  // Arrange GeoDivs by inset_pos
  std::map<std::string, std::vector<GeoDiv>> geo_divs_by_inset_pos;
  for (const InsetState &inset_state : inset_states_) {
    for (const auto &gd : inset_state.geo_divs()) {
      const std::string &id = gd.id();
      const std::string &inset_pos = csv_data.at(id).at("inset_pos");
      geo_divs_by_inset_pos[inset_pos].push_back(gd);
    }
  }

  // Create new InsetStates for each pos
  std::vector<InsetState> new_inset_states;
  for (auto &[inset_pos, geo_divs] : geo_divs_by_inset_pos) {
    InsetState new_inset_state(inset_pos, args_);
    for (auto &gd : geo_divs) {
      new_inset_state.push_back(std::move(gd));

      // Add target area, color, and label info to InsetState
      const std::string &id = gd.id();
      const auto &gd_info = csv_data.at(id);

      double target_area = std::stod(gd_info.at("area"));
      new_inset_state.insert_target_area(id, target_area);

      // Add color and label info, if present
      std::string color = gd_info.at("color");
      if (!color.empty())
        new_inset_state.insert_color(id, color);

      const std::string &label = gd_info.at("label");
      if (!label.empty())
        new_inset_state.insert_label(id, label);
    }
    new_inset_states.emplace_back(std::move(new_inset_state));
  }
  inset_states_ = std::move(new_inset_states);

  for (const auto &[id, data] : csv_data) {
    gd_to_inset_.emplace(id, data.at("inset_pos"));
  }
}

static bool is_point_as_separator(
  const std::map<std::string, std::map<std::string, std::string>> &csv_data)
{
  std::vector<std::string> area_strs;
  area_strs.reserve(csv_data.size());
  for (const auto &[_, data] : csv_data)
    area_strs.push_back(data.at("area"));
  return !StringToDecimalConverter::is_comma_as_separator(area_strs);
}

static void process_area_strs(
  std::map<std::string, std::map<std::string, std::string>> &csv_data)
{
  const bool uses_point_separator = is_point_as_separator(csv_data);
  for (auto &[_, data] : csv_data) {
    std::string &area_as_str = data.at("area");
    if (area_as_str.empty())
      area_as_str = "NA";
    area_as_str =
      StringToDecimalConverter::parse_str(area_as_str, uses_point_separator);
  }
}

void CartogramInfo::read_csv()
{
  ojson rows = load_csv_rows(args_.visual_file_name);
  const HeaderMap hm = headers_from(rows);

  if (hm.norm.size() < 2) {
    std::cerr
      << "ERROR: CSV with >= 2 columns (IDs, target areas) required. Some "
         "rows in your CSV may not have values for all columns"
      << std::endl;
    std::exit(17);
  }

  const std::string new_id_header = match_id_columns(args_.id_col);
  const int id_col = id_col_;  // index into hm.norm

  auto col_index = [&](const std::string &name) -> int {
    return hm.index_of_norm(name);
  };

  const int area_col = args_.area_col ? col_index(*args_.area_col) : 1;
  if (area_col < 0) {
    std::cerr
      << "ERROR: CSV with >= 2 columns (IDs, target areas) required. Some "
         "rows in your CSV may not have values for all columns"
      << std::endl;
    std::exit(17);
  }

  const int inset_col = col_index(args_.inset_col);  // default: "Inset"
  const int label_col = col_index(args_.label_col);  // default: "Label"
  const int color_col = extract_color_col_index(
    hm,
    args_.color_col);  // default: "Color" | "Colour"

  const std::string &id_key = hm.raw_at(static_cast<size_t>(id_col));
  const std::string &area_key = hm.raw_at(static_cast<size_t>(area_col));
  const std::string color_key = (color_col >= 0)
                                  ? hm.raw_at(static_cast<size_t>(color_col))
                                  : std::string();
  const std::string label_key = (label_col >= 0)
                                  ? hm.raw_at(static_cast<size_t>(label_col))
                                  : std::string();
  const std::string inset_key = (inset_col >= 0)
                                  ? hm.raw_at(static_cast<size_t>(inset_col))
                                  : std::string();

  std::map<std::string, std::map<std::string, std::string>> csv_data;

  for (const auto &r : rows.array_range()) {
    if (!r.contains(id_key) || !r.contains(area_key)) {
      std::cerr
        << "ERROR: CSV with >= 2 columns (IDs, target areas) required. Some "
           "rows in your CSV may not have values for all columns"
        << std::endl;
      std::exit(17);
    }

    std::string id = normalize_id_token(r.at(id_key).as_string());
    std::string area_as_str = r.at(area_key).as_string();
    check_validity_of_area_str(area_as_str);

    std::string color;
    if (!color_key.empty() && r.contains(color_key))
      color = r.at(color_key).as_string();

    std::string label;
    if (!label_key.empty() && r.contains(label_key))
      label = r.at(label_key).as_string();

    std::string inset_pos_as_str = "C";
    if (!inset_key.empty() && r.contains(inset_key))
      inset_pos_as_str = r.at(inset_key).as_string();

    const std::string inset_pos = process_inset_pos_str(inset_pos_as_str);
    check_validity_of_inset_pos(inset_pos, id);

    csv_data[id] = {
      {"area", area_as_str},
      {"color", color},
      {"label", label},
      {"inset_pos", inset_pos}};
  }

  update_id_header_info(new_id_header);
  check_validity_of_csv_ids(csv_data, initial_id_order_);
  process_area_strs(csv_data);
  relocate_geodivs_based_on_inset_pos(csv_data);
}
