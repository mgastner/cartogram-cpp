#ifndef CONFORMALITY_NODE_SET_HPP_
#define CONFORMALITY_NODE_SET_HPP_

#include "conformality/types.hpp"
#include <cstddef>
#include <vector>

namespace conformality
{

// Quadrature nodes that sample the current map T at fixed points of the
// original (equal-area) map. Stored as a structure of arrays; all vectors
// have the same length.
struct NodeSet {

  // Fixed original coordinates x_j
  std::vector<Vec2> x;

  // Quadrature weights: the original-frame area that node j represents
  std::vector<double> w;

  // Frozen Jacobian determinant g0_j = det DT0(x_j)
  std::vector<double> g0;

  // Current position y_j = T(x_j)
  std::vector<Vec2> y;

  // Current Jacobian DT(x_j)
  std::vector<Mat2> jac;

  [[nodiscard]] std::size_t size() const noexcept;

  // Throws std::invalid_argument if the vectors differ in length
  void check_consistent() const;
};

// Nodes at the centres of an nx-by-ny grid of equal cells covering
// [0, lx] x [0, ly], with T = identity (y = x, DT = I, g0 = 1) and each
// weight equal to the cell area. Throws std::invalid_argument for
// non-positive sizes.
[[nodiscard]] NodeSet make_cell_centre_nodes(
  double lx,
  double ly,
  std::size_t nx,
  std::size_t ny);

}  // namespace conformality

#endif  // CONFORMALITY_NODE_SET_HPP_
