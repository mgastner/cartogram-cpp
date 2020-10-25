#include "ft_real_2d.h"
#include <fftw3.h>
#include <iostream>

void FTReal2d::set_array_size(const unsigned int i, const unsigned int j)
{
  lx = i;
  ly = j;
  return;
}

void FTReal2d::ft_alloc()
{
  if (lx*ly <= 0) {
    std::cerr << "Invalid array dimensions in FTReal2dArray::ft_alloc()"
              << std::endl;
    _Exit(98915);
  }
  array = (double*) fftw_malloc(lx * ly * sizeof(double));
  return;
}

void FTReal2d::ft_free()
{
  if (array) {
    fftw_free(array);
    return;
  }
}

double *FTReal2d::get_array() const
{
  return array;
}

double &FTReal2d::operator ()(const unsigned int i, const unsigned int j)
{
  return array[i*ly + j];
}

const double FTReal2d::operator ()(const unsigned int i,
                                   const unsigned int j) const
{
  return array[i*ly + j];
}
