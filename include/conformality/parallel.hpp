#ifndef CONFORMALITY_PARALLEL_HPP_
#define CONFORMALITY_PARALLEL_HPP_

namespace conformality
{

// Whether the core was built with OpenMP (CMake option
// CONFORMALITY_OPENMP). Without it, every loop runs on one thread.
//
// The parallel loops never change results: the flow advances each point
// independently, and node sums are taken over fixed blocks of nodes whose
// partial sums are combined in block order. Results are therefore
// identical for every thread count, with or without OpenMP.
[[nodiscard]] bool openmp_enabled() noexcept;

// Number of threads used by the next parallel loop (1 without OpenMP)
[[nodiscard]] int max_threads() noexcept;

// Sets the number of threads for subsequent parallel loops. Has no effect
// without OpenMP or if n < 1.
void set_max_threads(int n) noexcept;

}  // namespace conformality

#endif  // CONFORMALITY_PARALLEL_HPP_
