#ifndef WRITE_TO_JSON_H_
#define WRITE_TO_JSON_H_

#include "cartogram_info.h"
#include "inset_state.h"
#include <nlohmann/json.hpp>
#include <string>
#include <iostream>
#include <boost/program_options.hpp>

nlohmann::json cgal_to_json(CartogramInfo *cart_info);
void write_to_json(nlohmann::json,
                   std::string,
                   std::string,
                   std::ostream&,
                   bool,
                   const boost::program_options::variables_map);
#endif
