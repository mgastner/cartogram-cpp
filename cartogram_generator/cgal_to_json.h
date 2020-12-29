#ifndef CGAL_TO_JSON_H_
#define CGAL_TO_JSON_H_

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "geo_div.h"

json cgal_to_json(std::vector<GeoDiv> container);

#endif
