#ifndef PROJECT_H_
#define PROJECT_H_

#include "cartogram_info.h"
#include "inset_state.h"

void project(InsetState*);
void project_with_triangulation(InsetState*);
void fill_graticule_diagonals(InsetState*);

#endif
