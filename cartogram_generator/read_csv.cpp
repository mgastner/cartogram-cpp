#include "csv.hpp"
#include <boost/program_options.hpp>
#include <iostream>

void read_csv(const boost::program_options::variables_map vm)
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
    std::cout << "ID is " << id << std::endl;

    // Get target area
    csv::CSVField area_field =
      vm.count("area") ? row[vm["area"].as<std::string>()] : row[1];
    double area = area_field.get<double>();
    std::cout << "Target area is " << area << std::endl;

    // Read color
    std::string color = "";
    if (vm.count("color")) {
      color = row[vm["color"].as<std::string>()].get();
    } else if (row.size() > 2) {
      color = row[2].get();
    }
    std::cout << "Color is " << color << std::endl;
  }
  return;
}
