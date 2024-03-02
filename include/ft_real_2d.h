#ifndef FT_REAL_2D_H_
#define FT_REAL_2D_H_

#include <fftw3.h>

class FTReal2d
{
private:
  double *array_ = nullptr;
  unsigned int lx_ = 0, ly_ = 0;  // Lattice dimensions
  fftw_plan plan_;

public:
  [[nodiscard]] double *as_1d_array() const;
  void set_array_size(unsigned int, unsigned int);
  void allocate(unsigned int, unsigned int);
  void free();
  void make_fftw_plan(const fftw_r2r_kind &, const fftw_r2r_kind &);
  void execute_fftw_plan();
  void destroy_fftw_plan();

  // Setter for array elements
  double &operator()(unsigned int, unsigned int);

  // Getter for array elements
  double operator()(unsigned int, unsigned int) const;
};

#endif
