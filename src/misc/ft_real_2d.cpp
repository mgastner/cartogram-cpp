#include "ft_real_2d.hpp"
#include <iostream>

double *FTReal2d::as_1d_array() const
{
  return array_;
}

void FTReal2d::set_array_size(const unsigned int i, const unsigned int j)
{
  lx_ = i;
  ly_ = j;
}

void FTReal2d::allocate(const unsigned int lx, const unsigned int ly)
{
  if (lx * ly <= 0) {
    std::cerr
      << "ERROR: Invalid array dimensions in FTReal2dArray::allocate_ft()"
      << std::endl;
    _Exit(98915);
  }
  lx_ = lx;
  ly_ = ly;
  array_ = static_cast<double *>(fftw_malloc(lx_ * ly_ * sizeof(double)));
}

void FTReal2d::free()
{
  fftw_free(array_);
}

void FTReal2d::make_fftw_plan(
  const fftw_r2r_kind &kind0,
  const fftw_r2r_kind &kind1)
{
  plan_ =
    fftw_plan_r2r_2d(lx_, ly_, array_, array_, kind0, kind1, FFTW_ESTIMATE);
}

void FTReal2d::execute_fftw_plan()
{
  fftw_execute(plan_);
}

void FTReal2d::destroy_fftw_plan()
{
  fftw_destroy_plan(plan_);
}

double &FTReal2d::operator()(const unsigned int i, const unsigned int j)
{
  return array_[i * ly_ + j];
}

double FTReal2d::operator()(const unsigned int i, const unsigned int j) const
{
  return array_[i * ly_ + j];
}
