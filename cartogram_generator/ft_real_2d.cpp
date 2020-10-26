#include "ft_real_2d.h"
#include <fftw3.h>
#include <iostream>

void FTReal2d::set_array_size(const unsigned int i, const unsigned int j)
{
  lx_ = i;
  ly_ = j;
  return;
}

void FTReal2d::allocate_ft()
{
  if (lx_*ly_ <= 0) {
    std::cerr << "Invalid array dimensions in FTReal2dArray::allocate_ft()"
              << std::endl;
    _Exit(98915);
  }
  array_ = (double*) fftw_malloc(lx_ * ly_ * sizeof(double));
  return;
}

void FTReal2d::free_ft()
{
  if (array_) {
    fftw_free(array_);
    return;
  }
}

const bool FTReal2d::is_allocated() const
{
  return !array_;
}

double *FTReal2d::array() const
{
  return array_;
}

double &FTReal2d::operator() (const unsigned int i, const unsigned int j)
{
  return array_[i*ly_ + j];
}

const double FTReal2d::operator() (const unsigned int i,
                                   const unsigned int j) const
{
  return array_[i*ly_ + j];
}
