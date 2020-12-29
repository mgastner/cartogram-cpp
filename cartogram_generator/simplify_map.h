#ifndef SIMPLIFY_MAP_H
#define SIMPLIFY_MAP_H

#include "cgal_typedef.h"
#include "map_state.h" 

#include <boost/program_options.hpp>
using namespace boost::program_options;
#include <boost/graph/adjacency_list.hpp>

void simplify_map(MapState *map_state);

#endif
