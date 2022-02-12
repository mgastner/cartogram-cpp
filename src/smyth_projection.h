#ifndef SMYTH_PROJECTION_H_
#define SMYTH_PROJECTION_H_

#include "inset_state.h"

double project_x_to_smyth(const double);
double project_y_to_smyth(const double);
void project_to_smyth_equal_surface(InsetState*);
void project_from_smyth_equal_surface(InsetState*);

#endif
