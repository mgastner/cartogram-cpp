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
  nlohmann::json j;

  // WHAT IF in_file isn't json???
  in_file >> j;

  return;
}
