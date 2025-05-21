#include "cartogram_info.hpp"
#include "csv.hpp"
#include "string_to_decimal_converter.hpp"

static int extract_color_col_index(
  const csv::CSVReader &reader,
  const std::string color_col_name)
{
  // Find index of column with color specifiers. If no color column header was
  // passed with the command-line flag --color, the header is assumed to be
  // "Color".
  int color_col = reader.index_of(color_col_name);

  // If the default "Color" header cannot be found, try again using the British
  // spelling "Colour"
  if (color_col == csv::CSV_NOT_FOUND && color_col_name == "Color") {
    color_col = reader.index_of("Colour");
  }

  return color_col;
}

static void check_validity_of_area_str(const std::string &area_as_str)
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
    std::exit(18);
  }

  if (
    !StringToDecimalConverter::is_str_NA(area_process_str) and
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

// Find the matching ID columns in both the CSV and GeoJSON file
// Returns the header name of the matching ID column in the GeoJSON file
std::string CartogramInfo::match_id_columns(
  const std::optional<std::string> &id_col)
{
  csv::CSVReader reader(args_.visual_file_name);
  std::string csv_id_header;

  // Copies what unique_properties_map_ does except it stores the values in a
  // set instead. Should find a more efficient way to do this.
  std::map<std::string, std::set<std::string>> geojson_properties_info;
  for (auto &[key, properties_vec] : unique_properties_map_) {
    std::set<std::string> properties_set(
      properties_vec.begin(),
      properties_vec.end());
    geojson_properties_info[key] = properties_set;
  }
  std::string matching_id_header;

  // If the user has specified a header as the ID column, check that one first.
  if (id_col) {
    std::set<std::string> csv_id_set;
    for (auto row = reader.begin(); row != reader.end(); row++)
      csv_id_set.insert((*row)[*id_col].get());

    // Check through each of the GEOJSON properties, see if any exactly match
    // the given CSV ID column data
    for (auto &[key, value_set] : geojson_properties_info) {
      if (csv_id_set == value_set) {
        matching_id_header = key;
        csv_id_header = *id_col;
        break;
      }
    }
    if (matching_id_header.empty())
      std::cerr << "Given ID header " << *id_col
                << " does not match with any GeoJSON properties. "
                   "Finding next best matching ID column..."
                << std::endl;
  }

  // If there is no user given ID header or the header does not match with any
  // GEOJSON properties, iterate through each of the CSV columns to find the
  // matching ID header.
  if (matching_id_header.empty()) {
    std::vector<std::string> column_headers = reader.get_col_names();

    for (std::string &header : column_headers) {
      // The begin() iterator for CSVReader seems to not be working correctly
      // (issue here: https://github.com/vincentlaucsb/csv-parser/issues/261).
      // As such, a new reader has to be declared in each loop in order to
      // properly iterate through the rows.
      csv::CSVReader loop_reader(args_.visual_file_name);
      std::set<std::string> data_set;
      for (auto row = loop_reader.begin(); row != loop_reader.end(); row++) {
        data_set.insert((*row)[header].get());
      }

      // If the set size is less than the number of rows then skip the column
      // as it cannot be the ID column.
      if (data_set.size() < loop_reader.n_rows())
        continue;

      for (auto &[key, value_set] : geojson_properties_info) {
        if (data_set == value_set) {
          matching_id_header = key;
          csv_id_header = header;
          break;
        }
      }

      if (!matching_id_header.empty())
        break;
    }
  }

  if (matching_id_header.empty()) {
    std::cerr << "ERROR: No valid matching ID header between GeoJSON "
                 "properties and CSV columns could be found."
              << std::endl;
    std::exit(16);
  }

  id_col_ = reader.index_of(csv_id_header);
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
    try {
      id = geojson_id_to_csv_id.at(id);
    } catch (const std::out_of_range &e) {
      std::cerr << "ERROR: Key '" << id
                << "' not found in geojson_id_to_csv_id. "
                << "Exception: " << e.what() << std::endl;
      // Re-throw, or return a default value
      throw;
    }
  }

  std::map<std::string, std::string> new_gd_to_inset;
  for (auto &[geojson_id, inset_pos] : gd_to_inset_) {
    try {
      const std::string csv_id = geojson_id_to_csv_id.at(geojson_id);
      new_gd_to_inset[csv_id] = inset_pos;
    } catch (const std::out_of_range &e) {
      std::cerr << "ERROR: Key '" << geojson_id
                << "' not found in geojson_id_to_csv_id. "
                << "Exception: " << e.what() << std::endl;
      // Re-throw, or return a default value
      throw;
    }
  }
  gd_to_inset_ = new_gd_to_inset;

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
  for (const auto &[id, data] : csv_data) {
    csv_ids.push_back(id);
  }

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
      std::cerr << "Warning: ID " << id << " in GeoJSON is not in CSV"
                << std::endl;
      csv_data[id] =
        {{"area", "NA"}, {"color", ""}, {"label", ""}, {"inset_pos", "C"}};
      // std::exit(22);
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
      std::map<std::string, std::string> gd_info;
      std::string inset_pos;
      try {
        gd_info = csv_data.at(id);
      } catch (const std::out_of_range &e) {
        std::cerr << "ERROR: Key '" << id << "' not found in csv_data. "
                  << "Exception: " << e.what() << std::endl;
        // Re-throw, or return a default value
        throw;
      }

      try {
        inset_pos = gd_info.at("inset_pos");
      } catch (const std::out_of_range &e) {
        std::cerr << "ERROR: Key '" << "inset_pos"
                  << "' not found in gd_info. "
                  << "Exception: " << e.what() << std::endl;
        // Re-throw, or return a default value
        throw;
      }
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
      std::map<std::string, std::string> gd_info;
      double target_area;
      std::string color;
      std::string label;

      try {
        gd_info = csv_data.at(id);
      } catch (const std::out_of_range &e) {
        std::cerr << "ERROR: Key '" << id << "' not found in csv_data. "
                  << "Exception: " << e.what() << std::endl;
        // Re-throw, or return a default value
        throw;
      }

      try {
        target_area = std::stod(gd_info.at("area"));
      } catch (const std::out_of_range &e) {
        std::cerr << "ERROR: Key '" << "area" << "' not found in gd_info. "
                  << "Exception: " << e.what() << std::endl;
        // Re-throw, or return a default value
        throw;
      }
      new_inset_state.insert_target_area(id, target_area);

      // Add color and label info, if present
      try {
        color = gd_info.at("color");
      } catch (const std::out_of_range &e) {
        std::cerr << "ERROR: Key '" << "color" << "' not found in gd_info. "
                  << "Exception: " << e.what() << std::endl;
        // Re-throw, or return a default value
        throw;
      }
      if (!color.empty()) {
        new_inset_state.insert_color(id, color);
      }

      try {
        label = gd_info.at("label");
      } catch (const std::out_of_range &e) {
        std::cerr << "ERROR: Key '" << "label" << "' not found in gd_info. "
                  << "Exception: " << e.what() << std::endl;
        // Re-throw, or return a default value
        throw;
      }
      if (!label.empty()) {
        new_inset_state.insert_label(id, label);
      }
    }
    new_inset_states.emplace_back(std::move(new_inset_state));
  }
  inset_states_ = std::move(new_inset_states);

  for (const auto &[id, data] : csv_data) {
    try {
      gd_to_inset_.emplace(id, data.at("inset_pos"));
    } catch (const std::out_of_range &e) {
      std::cerr << "ERROR: Key '" << "inset_pos" << "' not found in data. "
                << "Exception: " << e.what() << std::endl;
      // Re-throw, or return a default value
      throw;
    }
  }
}

static bool is_point_as_separator(
  const std::map<std::string, std::map<std::string, std::string>> &csv_data)
{
  std::vector<std::string> area_strs;
  for (const auto &[id, data] : csv_data) {
    try {
      area_strs.push_back(data.at("area"));
    } catch (const std::out_of_range &e) {
      std::cerr << "ERROR: Key '" << "area" << "' not found in data. "
                << "Exception: " << e.what() << std::endl;
      // Re-throw, or return a default value
      throw;
    }
  }

  if (StringToDecimalConverter::is_comma_as_separator(area_strs)) {
    return false;
  }

  return true;
}

static void process_area_strs(
  std::map<std::string, std::map<std::string, std::string>> &csv_data)
{
  const bool uses_point_separator = is_point_as_separator(csv_data);
  for (auto &[id, data] : csv_data) {
    try {
      std::string &area_as_str = data.at("area");
      if (area_as_str.empty()) {
        area_as_str = "NA";
      }
      area_as_str =
        StringToDecimalConverter::parse_str(area_as_str, uses_point_separator);
    } catch (const std::out_of_range &e) {
      std::cerr << "ERROR: Key '" << "area" << "' not found in data. "
                << "Exception: " << e.what() << std::endl;
      // Re-throw, or return a default value
      throw;
    }
  }
}

void CartogramInfo::read_csv()
{
  csv::CSVReader reader(args_.visual_file_name);

  const std::string new_id_header = match_id_columns(args_.id_col);
  const int id_col = id_col_;
  // Unless named through command-line argument,
  // 2nd column is assumed to be target areas
  const int area_col = args_.area_col ? reader.index_of(*args_.area_col) : 1;

  // Defaults set in parse_arguments.cpp
  const int inset_col = reader.index_of(args_.inset_col);  // default: "Inset"
  const int label_col = reader.index_of(args_.label_col);  // default: "Label"

  // default: "Color" | "Colour"
  const int color_col = extract_color_col_index(reader, args_.color_col);

  std::map<std::string, std::map<std::string, std::string>> csv_data;
  for (auto &row : reader) {
    if (row.size() < 2) {
      std::cerr
        << "ERROR: CSV with >= 2 columns (IDs, target areas) required. Some "
           "rows in your CSV may not have values for all columns"
        << std::endl;
      std::exit(17);
    }

    const std::string id = row[static_cast<size_t>(id_col)].get();
    const std::string area_as_str = row[static_cast<size_t>(area_col)].get();
    check_validity_of_area_str(area_as_str);

    const std::string color = (color_col != csv::CSV_NOT_FOUND)
                                ? row[static_cast<size_t>(color_col)].get()
                                : "";

    const std::string label = (label_col != csv::CSV_NOT_FOUND)
                                ? row[static_cast<size_t>(label_col)].get()
                                : "";

    const std::string inset_pos_as_str =
      (inset_col != csv::CSV_NOT_FOUND)
        ? row[static_cast<size_t>(inset_col)].get()
        : "C";

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
