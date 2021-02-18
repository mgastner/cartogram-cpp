#ifndef INTERPOL_H_
#define INTERPOL_H_

#include <boost/multi_array.hpp>

double interpolate_bilinearly(double x,
                              double y,
                              const boost::multi_array<double, 2> *grid,
                              char zero,
                              const unsigned int lx,
                              const unsigned int ly);

#endif
