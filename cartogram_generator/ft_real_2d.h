#ifndef FT_REAL_2D_H_
#define FT_REAL_2D_H_

#include <cstddef>

class FTReal2d {
  double *array_ = NULL;
  unsigned int lx_ = 0, ly_ = 0;    // Lattice dimensions
public:
  void set_array_size(const unsigned int, const unsigned int);
  void allocate_ft();
  void free_ft();
  const bool is_allocated() const;
  double *array() const;

  // Setter for array elements
  double &operator() (const unsigned int, const unsigned int);

  // Getter for array elements
  const double operator() (const unsigned int,
                           const unsigned int) const;
};

#endif
