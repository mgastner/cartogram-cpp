#ifndef INTERPOL_H_
#define INTERPOL_H_

#include <boost/multi_array.hpp>
#include <iostream>

double interpolate_bilinearly(
  double,
  double,
  boost::multi_array<double, 2> *,
  char,
  unsigned int,
  unsigned int);

double interpolate_bilinearly(
  double,
  double,
  std::function<double(unsigned int, unsigned int, char)> &,
  char,
  unsigned int,
  unsigned int);

#endif
