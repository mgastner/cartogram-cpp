#include "conformality/parallel.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace conformality
{

bool openmp_enabled() noexcept
{
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}

int max_threads() noexcept
{
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

void set_max_threads(const int n) noexcept
{
#ifdef _OPENMP
  if (n >= 1) {
    omp_set_num_threads(n);
  }
#else
  static_cast<void>(n);
#endif
}

}  // namespace conformality
