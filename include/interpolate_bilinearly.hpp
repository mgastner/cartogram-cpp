#ifndef INTERPOL_HPP_
#define INTERPOL_HPP_

#include <boost/multi_array.hpp>

double interpolate_bilinearly(
  double,
  double,
  const boost::multi_array<double, 2> &,
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

#endif  // INTERPOL_HPP_
