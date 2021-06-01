#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <map>

void print(std::map<std::string, std::vector<std::string> > properties_map)
{
  int i = 0;
  for (auto [key, value_vec] : properties_map)
  {
    i++;
    std::cout << i << ". " << key << ": ";
    for (std::string value : value_vec)
    {
      std::cout << value << ", ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
}

void geojson_to_csv(std::string geo_file_name)
{

  std::cout << "called geojson_to_csv.cpp" << std::endl;

  // Find all possible identifiers (key-value pairs where the
  // key is present in every feature and the value is unique for every feature)

  // Read input from a file
  std::ifstream in_file(geo_file_name);

  // Create nlohmann json object
  nlohmann::json j;
  in_file >> j;

  // {
  // NAME_1: [value, value, value],
  // GID_1: [value, value, value]
  // }

  // {GID_O: [BEL]}
  std::map<std::string, std::vector<std::string> > properties_map;
  for (auto feature : j["features"])
  {
    for (auto property_item : feature["properties"].items())
    {
      auto key = property_item.key();
      auto value = property_item.value();
      auto v = properties_map[key];
      bool value_not_inside = std::find(v.begin(), v.end(), value) == v.end();

      if (value != "" && value_not_inside)
      {
        properties_map[key].push_back(value);
      }
    }
  }

  std::map<std::string, std::vector<std::string> > viable_properties_map = properties_map;
  for (auto [key, value_vec] : properties_map)
  {
    if (value_vec.size() < j["features"].size())
    {
      viable_properties_map.erase(key);
    }
  }
  std::cout << std::endl;

  // Present user with all possible identifiers and a few examples
  std::cout << "Here are the available unique identifiers and their values. " << std::endl;
  print(viable_properties_map);
  std::cout << viable_properties_map.size() + 1 << ". All" << std::endl
            << std::endl;

  // Have the user choose which key they want to use as an identifier
  std::cout << "Please enter your number here: ";
  unsigned long chosen_number;
  std::cin >> chosen_number;
  std::cout << std::endl;

  // Declare chosen identifier
  std::map<std::string, std::vector<std::string>> chosen_identifiers;
  size_t i = 0;
  for (auto [key, value_vec] : viable_properties_map)
  {
    i++;
    if (chosen_number == i || chosen_number == viable_properties_map.size() + 1)
    {
      chosen_identifiers[key] = value_vec;
    }
  }
  print(chosen_identifiers);
}