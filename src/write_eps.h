#ifndef WRITE_EPS_H_
#define WRITE_EPS_H_

#include "inset_state.h"

void write_map_to_eps(const std::string, const bool, InsetState*);
void write_density_to_eps(const std::string, const double*, InsetState*);

#endif
