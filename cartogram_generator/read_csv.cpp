#include "map_state.h"
#include "csv.hpp"
#include <boost/program_options.hpp>
#include <iostream>

void read_csv(const boost::program_options::variables_map vm,
              CartogramInfo *cart_info)
{
  // Get name of CSV file from vm
  std::string csv_name;
  if (vm.count("visual_variable_file")) {
    csv_name = vm["visual_variable_file"].as<std::string>();
  } else {
    std::cerr << "ERROR: No CSV file given!" << std::endl;
    _Exit(16);
  }

  // Opening CSV Reader
  csv::CSVReader reader(csv_name);

  // Finding index of column with name "Inset"
  int inset_col = reader.index_of("Inset");

  // Reading CSV
  for (auto &row : reader) {
    if (row.size() < 2) {
      std::cerr << "ERROR: CSV with >= 2 columns (IDs, target areas) required"
                << std::endl;
      _Exit(17);
    }

    // Read ID of geographic division
    csv::CSVField id_field =
      vm.count("id") ? row[vm["id"].as<std::string>()] : row[0];
    std::string id = id_field.get();
    if (cart_info->ids_in_visual_variables_file().contains(id)) {
      std::cerr << "ERROR: ID "
                << id
                << " appears more than once in CSV"
                << std::endl;
      _Exit(301);
    }
    cart_info->insert_id_in_visual_variables_file(id);

    // Get target area
    csv::CSVField area_field =
      vm.count("area") ? row[vm["area"].as<std::string>()] : row[1];
    double area;
    if (!area_field.is_num()) {
      std::cout << "area_field" << area_field.get() << std::endl;
      if (area_field.get().compare("NA") == 0) {
        area = -1.0;  // Use negative area as sign of a missing value
      } else {
        std::cerr << "ERROR: Areas must be numeric or NA" << std::endl;
        _Exit(201);
      }
    } else {
      area = area_field.get<double>();
      if (area < 0.0) {
        std::cerr << "ERROR: negative area in CSV" << std::endl;
        _Exit(101);
      }
    }

    // Read color
    std::string color = "";
    if (vm.count("color")) {
      color = row[vm["color"].as<std::string>()].get();

    } else if (row.size() > 2) {
      color = row[2].get();
    }

    // Read inset
    std::string inset = "C"; // Assuming inset is C
    if (inset_col != csv::CSV_NOT_FOUND) {
      inset = row[inset_col].get();
    }

    InsetState inset_state(inset);
    inset_state.target_areas_insert(id, area);
    if (color != "") {
      inset_state.colors_insert(id, color);
    }

    cart_info->push_back(inset_state);

  }

  // Store header name of identifiers to read GeoJSON
  if (vm.count("id")) {
    cart_info->set_id_header(vm["id"].as<std::string>());
  } else {
    cart_info->set_id_header(reader.get_col_names()[0]);
  }
  return;
}
