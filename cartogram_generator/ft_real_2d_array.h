#ifndef FT_REAL_2D_ARRAY_H_
#define FT_REAL_2D_ARRAY_H_

#include <cstddef>

class FTReal2dArray {
  double *array = NULL;
public:
  double *get_array();
};

#endif
