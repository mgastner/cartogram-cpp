#ifndef WRITE_TO_JSON_H_
#define WRITE_TO_JSON_H_

#include "cartogram_info.h"
#include "inset_state.h"
#include <nlohmann/json.hpp>
#include <string>
#include <iostream>

nlohmann::json cgal_to_json(CartogramInfo *cart_info);
void write_to_json(nlohmann::json,
                   std::string,
                   std::string,
                   std::ostream&,
                   bool,
                   CartogramInfo *cart_info);
#endif
