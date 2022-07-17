#include "cartogram_info.h"
#include "csv.hpp"
#include <string>

void CartogramInfo::read_csv(const argparse::ArgumentParser &arguments)
{

  // Retrieve CSV Name.
  auto csv_name = arguments.get<std::string>("visual_variable_file");

  // Open CSV Reader.
  csv::CSVReader reader(csv_name);

  // Find index of column with IDs. If no ID column header was passed with the
  // command-line flag --id, the ID column is assumed to have index 0.
  auto is_id_header = arguments.present<std::string>("-D");
  int id_col = 0;
  if (is_id_header) {
    id_header_ = *is_id_header;
    id_col = reader.index_of(id_header_);
  } else {
    id_header_ = reader.get_col_names()[0];
  }

  // Find index of column with target areas. If no area column header was
  // passed with the command-line flag --area, the area column is assumed to
  // have index 1.
  int area_col = 1;
  if (auto area_header = arguments.present<std::string>("-A")) {
    std::cerr << "Area Header: " << *area_header << std::endl;
    area_col = reader.index_of(*area_header);
  }
  visual_variable_name_ = reader.get_col_names()[area_col];

  // Find index of column with inset specifiers. If no inset column header was
  // passed with the command-line flag --inset, the header is assumed to be
  // "Inset".
  auto inset_header = arguments.get<std::string>("-I");
  int inset_col = reader.index_of(inset_header);

  // Find index of column with color specifiers. If no color column header was
  // passed with the command-line flag --color, the header is assumed to be
  // "Color".
  auto color_header = arguments.get<std::string>("-C");
  int color_col = reader.index_of(color_header);

  // Default: "Label".
  auto label_header = arguments.get<std::string>("-L");
  int label_col = reader.index_of(label_header);

  // Read CSV
  std::set<std::string> inset_pos_set;
  for (auto &row : reader) {
    if (row.size() < 2) {
      std::cerr << "ERROR: CSV with >= 2 columns (IDs, target areas) required"
                << std::endl
                << "Some rows in your CSV may not have values for all columns"
                << std::endl;
      _Exit(17);
    }

    // Read ID of geographic division
    std::string id = row[id_col].get();
    if (ids_in_visual_variables_file_.contains(id)) {
      std::cerr << "ERROR: ID " << id << " appears more than once in CSV"
                << std::endl;
      _Exit(301);
    }
    ids_in_visual_variables_file_.insert(id);

    // Get target area
    csv::CSVField area_field = row[area_col];
    double area;
    if (area_field.is_num()) {
      area = area_field.get<double>();
      if (area < 0.0) {
        std::cerr << "ERROR: negative area in CSV" << std::endl;
        _Exit(101);
      }
    } else {
      std::string area_as_str = area_field.get();

      // With inspiration from:
      // https://stackoverflow.com/questions/2684491/remove-commas-from-string
      area_as_str.erase(
        std::remove(area_as_str.begin(), area_as_str.end(), ','),
        area_as_str.end());

      // Check if areas is missing or "NA"
      if (area_as_str.empty() || area_as_str == "NA") {
        area = -1.0;  // Use negative area as sign of a missing value
      }
      // With inspiration from:
      // https://stackoverflow.com/questions/4654636/how-to-determine-if-a-string-is-a-number-with-c
      else if (std::all_of(
                 area_as_str.begin(),
                 area_as_str.end(),
                 ::isdigit)) {
        area = std::stod(area_as_str);
      } else {
        std::cerr << "area_field: " << area_as_str << std::endl;
        std::cerr << "ERROR: Areas must be numeric or NA" << std::endl;
        _Exit(201);
      }
    }

    // Read color
    std::string color;
    if (color_col != csv::CSV_NOT_FOUND) {
      color = row[color_col].get();
    }

    // Read Label
    std::string label;
    if (label_col != csv::CSV_NOT_FOUND) {
      label = row[label_col].get();
    }

    // Read inset. Assume inset_pos is "C" if there is no inset column.
    std::string inset_pos = "C";
    if (inset_col != csv::CSV_NOT_FOUND) {
      inset_pos = row[inset_col].get();
      std::string inset_pos_original = inset_pos;

      // Set to "C" if inset position is blank
      if (inset_pos.empty()) {
        inset_pos = "C";
      }

      // Now we can process inputs like "center"/"left"/"right"
      inset_pos = std::toupper(inset_pos[0], std::locale());

      // Enable user to give inset position "U"/"D" for top and bottom inset
      if (inset_pos == "U") {
        inset_pos = "T";
      }
      if (inset_pos == "D") {
        inset_pos = "B";
      }

      // If unrecognized, set inset position to "C"
      std::unordered_set<std::string> permitted_pos{"C", "L", "R", "T", "B"};
      if (!permitted_pos.contains(inset_pos)) {
        std::cerr << "Unrecognized inset position : " << inset_pos_original
                  << " for Region: " << id << "\nSetting " << id
                  << "\'s inset position to Center (C)." << std::endl;
        inset_pos = "C";
      }
    }

    // Associate GeoDiv ID with inset position
    gd_to_inset_.insert(std::pair<std::string, std::string>(id, inset_pos));

    // Create inset_state for inset_pos unless it already exists
    if (!inset_pos_set.contains(inset_pos)) {
      InsetState inset_state(inset_pos);
      inset_states_.insert(
        std::pair<std::string, InsetState>(inset_pos, inset_state));
      inset_pos_set.insert(inset_pos);
    }

    // Insert target area and color
    InsetState *inset_state = &inset_states_.at(inset_pos);
    inset_state->insert_target_area(id, area);
    if (!color.empty()) {
      inset_state->insert_color(id, color);
    }
    if (!label.empty()) {
      inset_state->insert_label(id, label);
    }
  }
}
