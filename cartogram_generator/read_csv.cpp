#include <boost/program_options.hpp>
#include <iostream>

void read_csv(const boost::program_options::variables_map vm)
{
  // Get name of CSV EXIT_FAILURE
  std::string csv_name;
  if (vm.count("visual_variable_file")) {
    csv_name = vm["visual_variable_file"].as<std::string>();
  } else {
    std::cerr << "ERROR: No CSV file given!" << std::endl;
    _Exit(16);
  }

  return;
}
