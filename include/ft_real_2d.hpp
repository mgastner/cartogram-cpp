#ifndef FT_REAL_2D_HPP_
#define FT_REAL_2D_HPP_

#include <fftw3.h>

class FTReal2d
{
private:
  double *array_ = nullptr;
  unsigned int lx_ = 0, ly_ = 0;  // Lattice dimensions
  fftw_plan plan_;

public:
  [[nodiscard]] double *as_1d_array() const;
  void set_array_size(unsigned int i, unsigned int j);
  void allocate(unsigned int lx, unsigned int ly);
  void free();
  void make_fftw_plan(const fftw_r2r_kind &kind0, const fftw_r2r_kind &kind1);
  void execute_fftw_plan();
  void destroy_fftw_plan();

  // Setter for array elements
  double &operator()(unsigned int i, unsigned int j);

  // Getter for array elements
  double operator()(unsigned int i, unsigned int j) const;
};

#endif  // FT_REAL_2D_HPP_