#ifndef PROJECT_H_
#define PROJECT_H_

#include "map_state.h"

void project(MapState*);

void project_with_triangulation(MapState*);

void choose_diag(MapState*);

void choose_diag_2(MapState*);

void choose_diag_3(MapState*);

void choose_diag_4(MapState*);

void round_points(MapState*);

void point_search(MapState*, double, double, double, double);

#endif