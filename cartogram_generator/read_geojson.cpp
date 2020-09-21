#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>

void read_geojson(std::string geometry_file_name)
{
  std::ifstream in_file(geometry_file_name);
  if (!in_file) {
    throw std::system_error(errno,
                            std::system_category(),
                            "failed to open " + geometry_file_name);
  }

  using json = nlohmann::json;
  json j;
  try {
    in_file >> j;
  } catch (json::parse_error& e) {
    std::cerr << "ERROR: " << e.what() << '\n'
              << "exception id: " << e.id << '\n'
              << "byte position of error: " << e.byte << std::endl;
    _Exit(3);
  }

  return;
}
