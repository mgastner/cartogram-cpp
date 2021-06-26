// The class Real2dArray allows input of two-dimensional real-valued
// coordinates into a one-dimensional double array. The memory for a
// Real2dArray is allocated with fftw_malloc() instead of malloc() so that
// the object is optimized for Fourier transforms.

#ifndef REAL_2D_ARRAY_H_
#define REAL_2D_ARRAY_H_

#include <cstddef>
#include <fftw3.h>

class Real2dArray {
  double *array_ = NULL;
  unsigned int lx_ = 0, ly_ = 0;  // Lattice dimensions
  fftw_plan plan_;
public:
  double *as_1d_array();
  void allocate(const unsigned int, const unsigned int);
  void free();
  void make_fftw_plan();
  void execute_fftw_plan();
  void destroy_fftw_plan();

  // Setter for array elements
  double &operator() (const unsigned int, const unsigned int);

  // Getter for array elements
  double operator() (const unsigned int, const unsigned int) const;
};

#endif
