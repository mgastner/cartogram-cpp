#ifndef WRITE_TO_JSON_H_
#define WRITE_TO_JSON_H_

#include <string>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "cartogram_info.h"
#include "inset_state.h"

json cgal_to_json(InsetState*);

void write_to_json(json, std::string, std::string);

#endif
