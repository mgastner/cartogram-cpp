#ifndef CONFORMALITY_OBJECTIVE_HPP_
#define CONFORMALITY_OBJECTIVE_HPP_

#include "conformality/node_set.hpp"
#include "conformality/spectral_basis.hpp"
#include <vector>

namespace conformality
{

// Excess MIPS distortion of a 2x2 Jacobian with positive determinant,
//   m(J) = ||J||_F^2 / det J - 2 = K + 1/K - 2 >= 0,
// where K >= 1 is the ratio of the singular values. It vanishes exactly
// when J is a positive multiple of a rotation.
[[nodiscard]] double excess_mips(const Mat2 &jac) noexcept;

// Mean excess MIPS energy for alpha = 0 and unit importance weight,
//   E0 = sum_j w_j m(J_j) / sum_j w_j,
// with the actual det J_j in m. Returns +infinity if some det J_j <= 0
// (an inverted or degenerate node) or if the total weight is not positive.
// Uses compensated summation.
[[nodiscard]] double objective(const NodeSet &nodes);

// Modal derivatives F_mn = delta E0[T; v_mn], the first variation of E0
// along the velocity mode v_mn of `basis`, evaluated at the current nodes:
//
//   F_mn = (2 / Z) sum_j (w_j / det J_j) (J_j J_j^T) : Dv_mn(y_j),
//   Z = sum_j w_j.
//
// Because Dv_mn is traceless, this is the exact derivative of the discrete
// E0 (with actual det J) along the flow of v_mn. Returns a vector indexed
// like the basis coefficients. Throws std::domain_error if some
// det J_j <= 0 and std::invalid_argument if the total weight is not
// positive.
[[nodiscard]] std::vector<double> modal_derivatives(
  const NodeSet &nodes,
  const SpectralBasis &basis);

// Galerkin-projected H1 gradient and its descent direction
struct H1Descent {

  // Stream-function coefficients of the descent velocity, c_mn = -F_mn/G_mn
  std::vector<double> coeffs;

  // Squared H1 norm of the projected gradient, R = sum_mn F_mn^2 / G_mn.
  // The initial slope of E0 along the descent velocity is -R.
  double r{0.0};
};

// Computes c_mn = -F_mn / G_mn and R from modal derivatives f, which must
// be indexed like the coefficients of `basis`
[[nodiscard]] H1Descent h1_descent(
  const SpectralBasis &basis,
  const std::vector<double> &f);

}  // namespace conformality

#endif  // CONFORMALITY_OBJECTIVE_HPP_
