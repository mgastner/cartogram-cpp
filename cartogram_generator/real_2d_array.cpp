#include "real_2d_array.h"
#include <iostream>

double *Real2dArray::as_1d_array()
{
  return array_;
}

void Real2dArray::allocate(const unsigned int lx, const unsigned int ly)
{
  lx_ = lx;
  ly_ = ly;
  if (lx_*ly_ <= 0) {
    std::cerr << "Invalid array dimensions in Real2dArray::allocate()"
              << std::endl;
    _Exit(98915);
  }
  lx_ = lx;
  ly_ = ly;
  array_ = (double*) fftw_malloc(lx_ * ly_ * sizeof(double));
  return;
}

void Real2dArray::free()
{
  fftw_free(array_);
  return;
}

void Real2dArray::make_fftw_plan()
{
  plan_ = fftw_plan_r2r_2d(lx_, ly_,
                           array_, array_,
                           FFTW_RODFT01, FFTW_REDFT01, FFTW_ESTIMATE);
  return;
}

void Real2dArray::execute_fftw_plan()
{
  fftw_execute(plan_);
  return;
}

void Real2dArray::destroy_fftw_plan()
{
  fftw_destroy_plan(plan_);
  return;
}

// Setter for array elements
double &Real2dArray::operator() (const unsigned int i, const unsigned int j)
{
  return array_[i*ly_ + j];
}

// Getter for array elements
double Real2dArray::operator() (const unsigned int i,
                                const unsigned int j) const
{
  return array_[i*ly_ + j];
}
