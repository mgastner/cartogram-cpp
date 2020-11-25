#include "map_state.h"
#include "csv.hpp"
#include <boost/program_options.hpp>
#include <iostream>

void read_csv(const boost::program_options::variables_map vm,
              MapState *map_state)
{
  // Get name of CSV file from vm
  std::string csv_name;
  if (vm.count("visual_variable_file")) {
    csv_name = vm["visual_variable_file"].as<std::string>();
  } else {
    std::cerr << "ERROR: No CSV file given!" << std::endl;
    _Exit(16);
  }
  csv::CSVReader reader(csv_name);
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

    // Get target area
    csv::CSVField area_field =
      vm.count("area") ? row[vm["area"].as<std::string>()] : row[1];
    double area = area_field.get<double>();

    map_state->target_areas_insert(id, area);

    // Read color
    std::string color = "";
    if (vm.count("color")) {
      color = row[vm["color"].as<std::string>()].get();
      map_state->colors_insert(id, color);

    } else if (row.size() > 2) {
      color = row[2].get();
      map_state->colors_insert(id, color);
    }
  }

  // Store header name of identifiers to read GeoJSON
  if (vm.count("id")) {
    map_state->set_id_header(vm["id"].as<std::string>());
  } else {
    map_state->set_id_header(reader.get_col_names()[0]);
  }
  return;
}
