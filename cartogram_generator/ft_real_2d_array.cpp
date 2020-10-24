#include "ft_real_2d_array.h"
#include <fftw3.h>
#include <iostream>

void FTReal2dArray::set_array_size(const unsigned int i, const unsigned int j)
{
  lx = i;
  ly = j;
  return;
}

void FTReal2dArray::ft_alloc()
{
  if (lx*ly <= 0) {
    std::cerr << "Invalid array dimensions in FTReal2dArray::ft_alloc()"
              << std::endl;
    _Exit(98915);
  }
  array = (double*) fftw_malloc(lx * ly * sizeof(double));
  return;
}

void FTReal2dArray::ft_free()
{
  if (array) {
    fftw_free(array);
    return;
  }
}

double *FTReal2dArray::get_array()
{
  return array;
}

double &FTReal2dArray::operator ()(const unsigned int i, const unsigned int j)
{
  return array[i*ly + j];
}

double FTReal2dArray::operator ()(const unsigned int i,
                                  const unsigned int j) const
{
  return array[i*ly + j];
}
