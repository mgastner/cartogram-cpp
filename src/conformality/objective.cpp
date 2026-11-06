#include "conformality/objective.hpp"
#include <cassert>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace conformality
{

namespace
{

// Neumaier's compensated summation
class CompensatedSum
{
public:
  void add(const double value) noexcept
  {
    const double t = sum_ + value;
    if (std::abs(sum_) >= std::abs(value)) {
      compensation_ += (sum_ - t) + value;
    } else {
      compensation_ += (value - t) + sum_;
    }
    sum_ = t;
  }

  [[nodiscard]] double value() const noexcept
  {
    return sum_ + compensation_;
  }

private:
  double sum_{0.0};
  double compensation_{0.0};
};

double total_weight(const NodeSet &nodes)
{
  CompensatedSum z;
  for (const double w : nodes.w) {
    z.add(w);
  }
  return z.value();
}

}  // namespace

double excess_mips(const Mat2 &jac) noexcept
{
  return frobenius_sq(jac) / det(jac) - 2.0;
}

double objective(const NodeSet &nodes)
{
  nodes.check_consistent();
  const double z = total_weight(nodes);
  if (!(z > 0.0)) {
    return std::numeric_limits<double>::infinity();
  }
  CompensatedSum e;
  for (std::size_t j = 0; j < nodes.size(); ++j) {
    if (!(det(nodes.jac[j]) > 0.0)) {
      return std::numeric_limits<double>::infinity();
    }
    e.add(nodes.w[j] * excess_mips(nodes.jac[j]));
  }
  return e.value() / z;
}

std::vector<double> modal_derivatives(
  const NodeSet &nodes,
  const SpectralBasis &basis)
{
  nodes.check_consistent();
  const double z = total_weight(nodes);
  if (!(z > 0.0)) {
    throw std::invalid_argument(
      "modal_derivatives: total weight must be positive");
  }
  const std::size_t mx = basis.mx();
  const std::size_t my = basis.my();

  // With A = J J^T and S_m, C_m, T_n, U_n the sines and cosines of
  // kx_m y1 and ky_n y2,
  //   (A : Dv_mn) = (A11 - A22) kx_m ky_n C_m U_n
  //                 + A12 (kx_m^2 - ky_n^2) S_m T_n.
  // Accumulate the two node sums separately, then apply the wavenumbers.
  std::vector<double> cos_sum(mx * my, 0.0);
  std::vector<double> sin_sum(mx * my, 0.0);
  std::vector<double> sx(mx), cx(mx), sy(my), cy(my);
  for (std::size_t j = 0; j < nodes.size(); ++j) {
    const Mat2 &jac = nodes.jac[j];
    const double d = det(jac);
    if (!(d > 0.0)) {
      throw std::domain_error("modal_derivatives: det J must be positive");
    }
    const double a11 = jac.a11 * jac.a11 + jac.a12 * jac.a12;
    const double a22 = jac.a21 * jac.a21 + jac.a22 * jac.a22;
    const double a12 = jac.a11 * jac.a21 + jac.a12 * jac.a22;
    const double scale = 2.0 * nodes.w[j] / (z * d);
    const double a = scale * (a11 - a22);
    const double b = scale * a12;
    fill_harmonics(
      std::numbers::pi * nodes.y[j].x / basis.lx(),
      mx,
      sx.data(),
      cx.data());
    fill_harmonics(
      std::numbers::pi * nodes.y[j].y / basis.ly(),
      my,
      sy.data(),
      cy.data());
    for (std::size_t m = 0; m < mx; ++m) {
      const double am = a * cx[m];
      const double bm = b * sx[m];
      double *crow = cos_sum.data() + m * my;
      double *srow = sin_sum.data() + m * my;
      for (std::size_t n = 0; n < my; ++n) {
        crow[n] += am * cy[n];
        srow[n] += bm * sy[n];
      }
    }
  }

  std::vector<double> f(mx * my);
  for (std::size_t m = 1; m <= mx; ++m) {
    const double kx = basis.kx(m);
    for (std::size_t n = 1; n <= my; ++n) {
      const double ky = basis.ky(n);
      const std::size_t k = basis.index(m, n);
      f[k] = kx * ky * cos_sum[k] + (kx * kx - ky * ky) * sin_sum[k];
    }
  }
  return f;
}

H1Descent h1_descent(const SpectralBasis &basis, const std::vector<double> &f)
{
  if (f.size() != basis.n_modes()) {
    throw std::invalid_argument("h1_descent: wrong number of coefficients");
  }
  H1Descent result;
  result.coeffs.resize(f.size());
  CompensatedSum r;
  for (std::size_t m = 1; m <= basis.mx(); ++m) {
    for (std::size_t n = 1; n <= basis.my(); ++n) {
      const std::size_t k = basis.index(m, n);
      const double g = basis.h1_norm_sq(m, n);
      result.coeffs[k] = -f[k] / g;
      r.add(f[k] * f[k] / g);
    }
  }
  result.r = r.value();
  return result;
}

}  // namespace conformality
