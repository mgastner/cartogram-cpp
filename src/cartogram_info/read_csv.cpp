#include "cartogram_info.h"
#include "csv.hpp"
#include <string>

bool is_area_str_valid(const std::string &area_str)
{
  // Allow area_str being "NA"
  if (area_str == "NA") {
    return true;
  }

  // Only 0 to 9, '.', '-', and ',' are allowed
  for (const auto &c : area_str) {
    if (!(c >= '0' && c <= '9') && c != '.' && c != '-' && c != ',') {
      return false;
    }
  }

  // '-' can only be used as negative sign
  if (area_str.find('-') != std::string::npos) {
    if (area_str.find('-') != 0) {
      return false;
    }
  }
  return true;
}

void CartogramInfo::read_csv(const argparse::ArgumentParser &arguments)
{

  // Retrieve CSV Name
  const std::string csv_name =
    arguments.get<std::string>("visual_variable_file");

  // Open CSV Reader
  csv::CSVReader reader(csv_name);

  // Find index of column with IDs. If no ID column header was passed with the
  // command-line flag --id, the ID column is assumed to have index 0
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
  // have index 1
  int area_col = 1;
  if (auto area_header = arguments.present<std::string>("-A")) {
    std::cerr << "Area Header: " << *area_header << std::endl;
    area_col = reader.index_of(*area_header);
  }

  // Find index of column with inset specifiers. If no inset column header was
  // passed with the command-line flag --inset, the header is assumed to be
  // "Inset". This default value is set in parse_arguments.cpp.
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

    // Get target area as string
    csv::CSVField area_field = row[area_col];
    std::string area_as_str = area_field.get();

    // Parsed area string will be stored here
    double area;
    const char comma = ',';
    const char decimal = '.';
    const int area_str_size = area_as_str.length();

    /*
    Area String Parsing Algorithm:
    *) Check validity of string: only 0-9, ',' and '.' are allowed.

    *) Contains both decimal and comma (i.e: 123,456,789.123 or
    123.456.789,123):
      i) Comma appears first, decimal second (123,456,789.123):
      We consider it as US-convention, and remove the commas.
      123,456,789.123 -> 123456789.123.
      ii) Decimal appears first, comma second (123.456.789,123): We
      consider it Europe-convention, and we remove the decimals and replace the
      comma with decimal point (123.456.789,123 -> 123456789.123)
      iii) Number of comma and decimals both are more than 1
      (123,456,789.123.456): This format does not belong to any known
      convention, so exit with error.
      iv) Contains in random order (i.e: 123,456.789,123): Exit with error.

    *) Contains only commas (i.e: 123,456,789 or 123456,78):
      i) Only one comma is present and there are two digits after comma. We
      treat the comma as deicmal. (123456,78 -> 123456.78)
      ii) All other case, we assume the number does not contain any decimal
      values. So, we remove the commas and parse as usual. (123,456,789 ->
      123456789)

    *) Contains only decimals (i.e: 123.456.789 or 123456.789):
      i) Only one decimal present and there are three digits after that. Assume
      the decimal is European convention. (123456.789 -> 123456789)
      ii) Other cases, we assume the decimals were used to convey comma
      meaning. So, we remove the decimals places and parse the string.
      (123.456.789 -> 123456789)
    */

    if (!is_area_str_valid(area_as_str)) {
      std::cerr << "ERROR: Invalid area string: " << area_as_str << std::endl;
      std::cerr << "Area string must only contain 0-9, '.', '-' and ','."
                << std::endl;
      _Exit(18);
    }

    // Both commas and decimals are present
    if (
      (area_as_str.find(comma) != std::string::npos) &&
      (area_as_str.find(decimal) != std::string::npos)) {

      // Number of comma or decimals both more than 1, exit
      // handles 123,456,789.123.456 case
      int comma_count =
        std::count(area_as_str.begin(), area_as_str.end(), comma);
      int decimal_count =
        std::count(area_as_str.begin(), area_as_str.end(), decimal);
      if (comma_count > 1 and decimal_count > 1) {
        std::cerr << "ERROR: Invalid area string: " << area_as_str
                  << std::endl;
        _Exit(19);
      }

      // TODO: Add logic to handle cases like 123,456.789,123
      // and 123.456,789.123

      // if commas come first, remove all commas and keep the decimal
      if (area_as_str.find(comma) < area_as_str.find(decimal)) {
        area_as_str.erase(
          std::remove(area_as_str.begin(), area_as_str.end(), comma),
          area_as_str.end());
      } else {  // decimal comes first

        // remove all decimals and replace commas with decimal
        area_as_str.erase(
          std::remove(area_as_str.begin(), area_as_str.end(), decimal),
          area_as_str.end());
        area_as_str[area_as_str.find(comma)] = '.';
      }
      area = std::stod(area_as_str);
    } else if (area_as_str.find(comma) != std::string::npos) {  // only commas
      int comma_count =
        std::count(area_as_str.begin(), area_as_str.end(), comma);

      // if there is only 1 comma, and there are only 2 numbers after the
      // comma, convert to decimal: 12,56 -> 12.56
      if (
        (comma_count == 1) &&
        (area_str_size - 1 - area_as_str.find(comma) == 2)) {
        area_as_str[area_as_str.find(comma)] = '.';
      } else {

        // treat as US-convention, remove commas
        area_as_str.erase(
          std::remove(area_as_str.begin(), area_as_str.end(), comma),
          area_as_str.end());
      }
      area = std::stod(area_as_str);
    } else if (area_as_str.find(decimal) != std::string::npos) {  // only
                                                                  // decimals
      int decimal_count =
        std::count(area_as_str.begin(), area_as_str.end(), decimal);

      // If there is more than 1 decimal, remove all decimals, as we treat them
      // as commas
      if (decimal_count > 1) {
        area_as_str.erase(
          std::remove(area_as_str.begin(), area_as_str.end(), decimal),
          area_as_str.end());
      } else if (
        (decimal_count == 1) &&
        (area_str_size - 1 - area_as_str.find(decimal) == 3)) {

        // if there is only 1 decimal, and there are only 3 numbers after the
        // decimal, treat it as European-convention, pretend it is a comma.
        // 12456.789 -> 12456789
        area_as_str.erase(
          std::remove(area_as_str.begin(), area_as_str.end(), decimal),
          area_as_str.end());
      }
      area = std::stod(area_as_str);
    } else if (area_as_str.empty() || area_as_str.compare("NA") == 0) {
      area = -1.0;  // Use negative area as sign of a missing value
    } else {
      area = std::stod(area_as_str);
    }

    if (
      area < 0.0 and
      (area_as_str.compare("NA") != 0 and !area_as_str.empty())) {
      std::cerr << "ERROR: Negative area in CSV" << std::endl;
      _Exit(101);
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
