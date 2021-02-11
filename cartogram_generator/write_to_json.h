#ifndef WRITE_TO_JSON_H_
#define WRITE_TO_JSON_H_

#include <string>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "map_state.h"

json cgal_to_json(MapState *map_state);

void write_to_json(json container, std::string geo_file_name, std::string new_geo_fn);

#endif