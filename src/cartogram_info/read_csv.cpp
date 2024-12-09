#include "cartogram_info.hpp"
#include "csv.hpp"
#include "string_to_decimal_converter.hpp"

csv::CSVReader load_csv(const argparse::ArgumentParser &arguments)
{
  // Retrieve CSV Name
  const std::string csv_filename =
    arguments.get<std::string>("visual_variable_file");

  // Open CSV Reader
  csv::CSVReader reader(csv_filename);

  return reader;
}

int extract_id_header_col_index(
  const csv::CSVReader &reader,
  const argparse::ArgumentParser &arguments)
{
  // Find index of column with IDs. If no ID column header was passed with the
  // command-line flag --id, the ID column is assumed to have index 0
  int id_col = 0;
  if (auto is_id_header = arguments.present<std::string>("-D")) {
    const std::string id_header = *is_id_header;
    id_col = reader.index_of(id_header);
  }

  return id_col;
}

int extract_area_col_index(
  const csv::CSVReader &reader,
  const argparse::ArgumentParser &arguments)
{
  // Find index of column with target areas. If no area column header was
  // passed with the command-line flag --area, the area column is assumed to
  // have index 1
  int area_col = 1;
  if (auto area_header = arguments.present<std::string>("-A")) {
    area_col = reader.index_of(*area_header);
  }

  return area_col;
}

int extract_inset_col_index(
  const csv::CSVReader &reader,
  const argparse::ArgumentParser &arguments)
{
  // Find index of column with inset specifiers. If no inset column header was
  // passed with the command-line flag --inset, the header is assumed to be
  // "Inset". This default value is set in parse_arguments.cpp.
  auto inset_header = arguments.get<std::string>("-I");
  int inset_col = reader.index_of(inset_header);

  return inset_col;
}

int extract_color_col_index(
  const csv::CSVReader &reader,
  const argparse::ArgumentParser &arguments)
{
  // Find index of column with color specifiers. If no color column header was
  // passed with the command-line flag --color, the header is assumed to be
  // "Color".
  auto color_header = arguments.get<std::string>("-C");
  int color_col = reader.index_of(color_header);

  return color_col;
}

int extract_label_col_index(
  const csv::CSVReader &reader,
  const argparse::ArgumentParser &arguments)
{
  // Default: "Label".
  auto label_header = arguments.get<std::string>("-L");
  int label_col = reader.index_of(label_header);

  return label_col;
}

void check_validity_of_area_str(const std::string &area_as_str)
{
  std::string area_process_str = area_as_str;

  if (area_process_str.empty()) {
    area_process_str = "NA";
  }

  if (!StringToDecimalConverter::is_str_valid_characters(area_process_str)) {
    std::cerr
      << "ERROR: Invalid area string: " << area_process_str
      << ". Area string must only contain 0-9, '.', '-' and ',' or 'NA'."
      << std::endl;
    _Exit(18);
  }

  if (
    !StringToDecimalConverter::is_str_NA(area_process_str) and
    !StringToDecimalConverter::is_str_correct_format(area_process_str)) {
    std::cerr << "ERROR: Invalid area string format: " << area_process_str
              << std::endl;
    _Exit(19);
  }

  if (area_process_str.front() == '-') {
    std::cerr << "ERROR: Negative area in CSV" << std::endl;
    _Exit(101);
  }
}

std::string process_inset_pos_str(const std::string &inset_pos_as_str)
{
  std::string inset_pos = inset_pos_as_str;

  if (inset_pos.empty()) {
    inset_pos = "C";
  }

  inset_pos = std::toupper(inset_pos[0], std::locale());

  if (inset_pos == "U") {
    inset_pos = "T";
  }
  if (inset_pos == "D") {
    inset_pos = "B";
  }
  return inset_pos;
}

void check_validity_of_inset_pos(
  const std::string &inset_pos,
  const std::string &id)
{
  std::unordered_set<std::string> permitted_pos{"C", "L", "R", "T", "B"};
  if (!permitted_pos.contains(inset_pos)) {
    std::cerr << "Unrecognized inset position : " << inset_pos
              << " for Region: " << id << "\nSetting " << id
              << "\'s inset position to Center (C)." << std::endl;
    _Exit(20);
  }
}

void CartogramInfo::update_id_header_info(const std::string &csv_id_header)
{
  if (
    std::find(
      unique_properties_.begin(),
      unique_properties_.end(),
      csv_id_header) == unique_properties_.end()) {
    std::cerr
      << "ERROR: ID header not found in GeoJSON properties or is not unique."
      << std::endl;
    std::cerr << "Provided ID header: " << csv_id_header << std::endl;
    std::cerr << "Unique properties in GeoJSON: ";
    for (const auto &prop : unique_properties_) {
      std::cerr << prop << ", ";
    }
    _Exit(16);
  }

  // The automically chosen ID header from read GeoJSON is the same as the ID
  // header from read CSV
  if (csv_id_header == id_header_) {
    return;
  }

  std::map<std::string, std::string> geojson_id_to_csv_id;
  for (auto &[geojson_id, properties] : properties_map_) {
    const std::string csv_id = properties.at(csv_id_header);
    geojson_id_to_csv_id[geojson_id] = csv_id;
  }

  for (auto &id : initial_id_order_) {
    id = geojson_id_to_csv_id.at(id);
  }

  std::map<std::string, std::string> new_gd_to_inset;
  for (auto &[geojson_id, inset_pos] : gd_to_inset_) {
    const std::string csv_id = geojson_id_to_csv_id.at(geojson_id);
    new_gd_to_inset[csv_id] = inset_pos;
  }
  gd_to_inset_ = new_gd_to_inset;

  for (auto &[inset_pos, inset_state] : inset_states_) {
    inset_state.update_gd_ids(geojson_id_to_csv_id);
  }

  std::map<std::string, std::map<std::string, std::string>> new_properties_map;
  for (auto &[geojson_id, properties] : properties_map_) {
    const std::string csv_id = geojson_id_to_csv_id.at(geojson_id);
    new_properties_map[csv_id] = properties;
  }
  properties_map_ = new_properties_map;

  id_header_ = csv_id_header;
}

void check_validity_of_csv_ids(
  const std::map<std::string, std::map<std::string, std::string>> &csv_data,
  const std::vector<std::string> &initial_id_order)
{
  std::vector<std::string> csv_ids;
  for (const auto &[id, data] : csv_data) {
    csv_ids.push_back(id);
  }

  for (const auto &id : csv_ids) {
    if (
      std::find(initial_id_order.begin(), initial_id_order.end(), id) ==
      initial_id_order.end()) {
      std::cerr << "ERROR: ID " << id << " in CSV is not in GeoJSON"
                << std::endl;
      _Exit(21);
    }
  }

  for (const auto &id : initial_id_order) {
    if (std::find(csv_ids.begin(), csv_ids.end(), id) == csv_ids.end()) {
      std::cerr << "ERROR: ID " << id << " in GeoJSON is not in CSV"
                << std::endl;
      _Exit(22);
    }
  }
}

void CartogramInfo::relocate_geodivs_based_on_inset_pos(
  const std::map<std::string, std::map<std::string, std::string>> &csv_data)
{
  std::map<std::string, InsetState> new_inset_states;
  for (const auto &[inset_pos, inset_state] : inset_states_) {
    for (const auto &gd : inset_state.geo_divs()) {
      const std::string id = gd.id();
      const std::string inset_pos = csv_data.at(id).at("inset_pos");
      if (!new_inset_states.contains(inset_pos)) {
        InsetState new_inset_state(inset_pos);
        new_inset_states.insert(
          std::pair<std::string, InsetState>(inset_pos, new_inset_state));
      }
      InsetState *new_inset_state = &new_inset_states.at(inset_pos);
      new_inset_state->push_back(gd);
      new_inset_state->insert_target_area(
        id,
        std::stod(csv_data.at(id).at("area")));
      if (!csv_data.at(id).at("color").empty()) {
        std::string color = csv_data.at(id).at("color");
        new_inset_state->insert_color(id, color);
      }
      if (!csv_data.at(id).at("label").empty()) {
        std::string label = csv_data.at(id).at("label");
        new_inset_state->insert_label(id, label);
      }
    }
  }
  inset_states_ = new_inset_states;

  for (const auto &[id, data] : csv_data) {
    gd_to_inset_.insert(
      std::pair<std::string, std::string>(id, data.at("inset_pos")));
  }
}

bool is_point_as_separator(
  const std::map<std::string, std::map<std::string, std::string>> &csv_data)
{
  std::vector<std::string> area_strs;
  for (const auto &[id, data] : csv_data) {
    area_strs.push_back(data.at("area"));
  }

  if (StringToDecimalConverter::is_comma_as_separator(area_strs)) {
    return false;
  }

  return true;
}

void process_area_strs(
  std::map<std::string, std::map<std::string, std::string>> &csv_data)
{
  const bool uses_point_separator = is_point_as_separator(csv_data);
  for (auto &[id, data] : csv_data) {
    std::string &area_as_str = data.at("area");
    if (area_as_str.empty()) {
      area_as_str = "NA";
    }
    area_as_str =
      StringToDecimalConverter::parse_str(area_as_str, uses_point_separator);
  }
}

void CartogramInfo::read_csv(const argparse::ArgumentParser &arguments)
{
  csv::CSVReader reader = load_csv(arguments);

  const int id_col = extract_id_header_col_index(reader, arguments);
  const int area_col = extract_area_col_index(reader, arguments);
  const int inset_col = extract_inset_col_index(reader, arguments);
  const int color_col = extract_color_col_index(reader, arguments);
  const int label_col = extract_label_col_index(reader, arguments);

  std::map<std::string, std::map<std::string, std::string>> csv_data;
  for (auto &row : reader) {
    if (row.size() < 2) {
      std::cerr
        << "ERROR: CSV with >= 2 columns (IDs, target areas) required. Some "
           "rows in your CSV may not have values for all columns"
        << std::endl;
      _Exit(17);
    }

    const std::string id = row[id_col].get();
    const std::string area_as_str = row[area_col].get();
    check_validity_of_area_str(area_as_str);

    const std::string color =
      (color_col != csv::CSV_NOT_FOUND) ? row[color_col].get() : "";

    const std::string label =
      (label_col != csv::CSV_NOT_FOUND) ? row[label_col].get() : "";

    const std::string inset_pos_as_str =
      (inset_col != csv::CSV_NOT_FOUND) ? row[inset_col].get() : "C";

    const std::string inset_pos = process_inset_pos_str(inset_pos_as_str);
    check_validity_of_inset_pos(inset_pos, id);

    csv_data[id] = {
      {"area", area_as_str},
      {"color", color},
      {"label", label},
      {"inset_pos", inset_pos}};
  }

  process_area_strs(csv_data);

  const std::string id_header = reader.get_col_names()[id_col];
  update_id_header_info(id_header);
  check_validity_of_csv_ids(csv_data, initial_id_order_);
  relocate_geodivs_based_on_inset_pos(csv_data);
}
