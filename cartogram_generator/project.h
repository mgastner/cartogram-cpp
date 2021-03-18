#ifndef PROJECT_H_
#define PROJECT_H_

#include "map_state.h"

void project(MapState*);

std::vector<int> find_graticule_point(const int,
                                      const int,
                                      const int,
                                      const int);

void project_graticule_centroids(MapState*);

void project_with_triangulation(MapState*);

void choose_diag(MapState*);

#endif