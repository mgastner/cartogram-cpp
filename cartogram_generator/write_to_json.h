#ifndef WRITE_TO_JSON_H_
#define WRITE_TO_JSON_H_

#include <string>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "cartogram_info.h"
#include "inset_state.h"
#include <iostream>

json cgal_to_json(InsetState*);
void write_to_json(json, std::string, std::string, CGAL::Bbox_2);
json cgal_to_json_all_insets(CartogramInfo *cart_info);
void write_to_json_all_insets(json, std::string, std::ostream&);
// void write_to_json_all_frames(json, std::string,
//                               std::string,
//                               std::map <std::string, CGAL::Bbox_2>);

#endif
