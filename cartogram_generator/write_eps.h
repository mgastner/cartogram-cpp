#ifndef WRITE_EPS_H_
#define WRITE_EPS_H_

void write_map_to_eps(std::string, InsetState*);
void write_density_to_eps(std::string, Real2dArray, InsetState*);

void write_velocity_to_eps(std::string,
                           boost::multi_array<double, 2>*,
                           InsetState*);
#endif
