#ifndef WRITE_TO_JSON_H_
#define WRITE_TO_JSON_H_

#include <string>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

void write_to_json(json oldJ, json container, std::string geo_file_name);

#endif
