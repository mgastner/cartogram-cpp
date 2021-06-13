#ifndef PROJECT_H_
#define PROJECT_H_

#include "cartogram_info.h"
#include "inset_state.h"

void project(InsetState*);

void project_with_triangulation(InsetState*);

void choose_diag(InsetState*);

void choose_diag_2(InsetState*);

void choose_diag_3(InsetState*);

void choose_diag_4(InsetState*);

void round_points(InsetState*);

void point_search(InsetState*, double, double, double, double);

#endif