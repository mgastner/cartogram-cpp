// TODO: Import optional color ID.
#include <iostream>
#include "read_csv.h"
#include <boost/program_options.hpp>
#include "csv.hpp"

using namespace boost::program_options;
using namespace csv;

void read_csv(const variables_map vm) {
  std::cout << "In read_csv()" << std::endl;

  std::string csv_name;

  if (vm.count("visual_variables")) {
    csv_name = vm["visual_variables"].as<std::string>();
  }
  else {
    std::cout << "No CSV file given! \n\n";
    // exit();
  }

  CSVReader reader(csv_name);

  // to store the data, to change to map next
  std::vector<std::string> region_names;
  std::vector<double> target_areas;
  std::vector<std::string> colors;

  for (auto& row: reader) {

    if (row.size() < 2) {
      break;
      std::cout << "A CSV with at least 2 columns (Identifiers, Target Areas) is required";
      // exit();
    }

    // check whether flags are present in variable map
    if (vm.count("id")) {
      region_names.push_back(row[vm["id"].as<std::string>()].get());
    }
    else {
      region_names.push_back(row[0].get());
    }

    if (vm.count("area")) {
      target_areas.push_back(row[vm["area"].as<std::string>()].get<double>());
    }
    else {
      target_areas.push_back(row[1].get<double>());
    }

    if (vm.count("color")) {
      colors.push_back(row[vm["color"].as<std::string>()].get());
    }
    else if (row.size() > 2){
      colors.push_back(row[2].get());
    }

  }

  // printing data for now
  for (size_t i = 0; i < region_names.size(); i++) {
    std::cout << "Name: " << region_names[i] << '\n';
    std::cout << "Target Area: " << target_areas[i] << '\n';
    std::cout << "Color: " << colors[i] << "\n \n";
  }


  return;
}
