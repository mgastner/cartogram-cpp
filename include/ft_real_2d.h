#ifndef FT_REAL_2D_H_
#define FT_REAL_2D_H_

#include <fftw3.h>
#include <cstddef>

class FTReal2d {
private:
  double *array_ = NULL;
  unsigned int lx_ = 0, ly_ = 0;  // Lattice dimensions
  fftw_plan plan_;

public:
  double *as_1d_array() const;
  void set_array_size(const unsigned int, const unsigned int);
  void allocate(const unsigned int, const unsigned int);
  void free();
  void make_fftw_plan(fftw_r2r_kind, fftw_r2r_kind);
  void execute_fftw_plan();
  void destroy_fftw_plan();

  // Setter for array elements
  double &operator() (const unsigned int, const unsigned int);

  // Getter for array elements
  double operator() (const unsigned int, const unsigned int) const;
};

#endif
