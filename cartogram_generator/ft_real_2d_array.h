#ifndef FT_REAL_2D_ARRAY_H_
#define FT_REAL_2D_ARRAY_H_

#include <cstddef>

class FTReal2dArray {
  double *array = NULL;
  unsigned int lx = 0, ly = 0;  // Lattice dimensions
public:
  void set_array_size(const unsigned int, const unsigned int);
  void ft_alloc();
  void ft_free();
  double *get_array();
  double &operator ()(const unsigned int, const unsigned int);  // Setter
  double operator ()(const unsigned int, const unsigned int) const;  // Getter
};

#endif
