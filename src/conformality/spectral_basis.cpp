#include "conformality/spectral_basis.hpp"
#include <cassert>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace conformality
{

void fill_harmonics(
  const double theta,
  const std::size_t count,
  double *sin_out,
  double *cos_out) noexcept
{
  if (count == 0) {
    return;
  }
  const double s1 = std::sin(theta);
  const double c1 = std::cos(theta);
  sin_out[0] = s1;
  cos_out[0] = c1;
  for (std::size_t j = 1; j < count; ++j) {

    // sin((j + 1) theta) and cos((j + 1) theta) from the values at j theta
    sin_out[j] = sin_out[j - 1] * c1 + cos_out[j - 1] * s1;
    cos_out[j] = cos_out[j - 1] * c1 - sin_out[j - 1] * s1;
  }
}

SpectralBasis::SpectralBasis(
  const double lx,
  const double ly,
  const std::size_t mx,
  const std::size_t my)
    : lx_(lx), ly_(ly), mx_(mx), my_(my)
{
  if (!(lx > 0.0) || !(ly > 0.0)) {
    throw std::invalid_argument("SpectralBasis: lx and ly must be positive");
  }
  if (mx == 0 || my == 0) {
    throw std::invalid_argument("SpectralBasis: mx and my must be positive");
  }
  kx_.resize(mx_);
  kx_sq_.resize(mx_);
  for (std::size_t m = 1; m <= mx_; ++m) {
    const double k = static_cast<double>(m) * std::numbers::pi / lx_;
    kx_[m - 1] = k;
    kx_sq_[m - 1] = k * k;
  }
  ky_.resize(my_);
  ky_sq_.resize(my_);
  for (std::size_t n = 1; n <= my_; ++n) {
    const double k = static_cast<double>(n) * std::numbers::pi / ly_;
    ky_[n - 1] = k;
    ky_sq_[n - 1] = k * k;
  }
}

double SpectralBasis::lx() const noexcept
{
  return lx_;
}

double SpectralBasis::ly() const noexcept
{
  return ly_;
}

std::size_t SpectralBasis::mx() const noexcept
{
  return mx_;
}

std::size_t SpectralBasis::my() const noexcept
{
  return my_;
}

std::size_t SpectralBasis::n_modes() const noexcept
{
  return mx_ * my_;
}

std::size_t SpectralBasis::index(const std::size_t m, const std::size_t n)
  const noexcept
{
  assert(m >= 1 && m <= mx_ && n >= 1 && n <= my_);
  return (m - 1) * my_ + (n - 1);
}

double SpectralBasis::kx(const std::size_t m) const noexcept
{
  assert(m >= 1 && m <= mx_);
  return kx_[m - 1];
}

double SpectralBasis::ky(const std::size_t n) const noexcept
{
  assert(n >= 1 && n <= my_);
  return ky_[n - 1];
}

const std::vector<double> &SpectralBasis::kx_all() const noexcept
{
  return kx_;
}

const std::vector<double> &SpectralBasis::kx_sq_all() const noexcept
{
  return kx_sq_;
}

const std::vector<double> &SpectralBasis::ky_all() const noexcept
{
  return ky_;
}

const std::vector<double> &SpectralBasis::ky_sq_all() const noexcept
{
  return ky_sq_;
}

double SpectralBasis::h1_norm_sq(const std::size_t m, const std::size_t n)
  const noexcept
{
  const double k_sq = kx_sq_[m - 1] + ky_sq_[n - 1];
  return 0.25 * lx_ * ly_ * k_sq * k_sq;
}

SpectralBasis::Field SpectralBasis::evaluate(
  const std::vector<double> &c,
  const Vec2 y,
  Workspace &ws) const
{
  assert(c.size() == n_modes());
  ws.sin_x.resize(mx_);
  ws.cos_x.resize(mx_);
  ws.sin_y.resize(my_);
  ws.cos_y.resize(my_);
  fill_harmonics(
    std::numbers::pi * y.x / lx_,
    mx_,
    ws.sin_x.data(),
    ws.cos_x.data());
  fill_harmonics(
    std::numbers::pi * y.y / ly_,
    my_,
    ws.sin_y.data(),
    ws.cos_y.data());

  // With S_m = sin(kx_m y1), C_m = cos(kx_m y1), T_n = sin(ky_n y2), and
  // U_n = cos(ky_n y2), the inner sums over n are
  //   p_m = sum_n c_mn T_n,  q_m = sum_n c_mn ky_n U_n,
  //   r_m = sum_n c_mn ky_n^2 T_n,
  // and then
  //   s = sum_m S_m p_m,             v1 = sum_m S_m q_m,
  //   v2 = -sum_m kx_m C_m p_m,      dv11 = sum_m kx_m C_m q_m,
  //   dv12 = -sum_m S_m r_m,         dv21 = sum_m kx_m^2 S_m p_m,
  //   dv22 = -dv11.
  Field f;
  for (std::size_t i = 0; i < mx_; ++i) {
    const double *row = c.data() + i * my_;
    double p = 0.0;
    double q = 0.0;
    double r = 0.0;
    for (std::size_t j = 0; j < my_; ++j) {
      const double ct = row[j] * ws.sin_y[j];
      p += ct;
      q += row[j] * ky_[j] * ws.cos_y[j];
      r += ct * ky_sq_[j];
    }
    const double sm = ws.sin_x[i];
    const double kcm = kx_[i] * ws.cos_x[i];
    f.s += sm * p;
    f.v.x += sm * q;
    f.v.y -= kcm * p;
    f.dv.a11 += kcm * q;
    f.dv.a12 -= sm * r;
    f.dv.a21 += kx_sq_[i] * sm * p;
  }
  f.dv.a22 = -f.dv.a11;
  return f;
}

SpectralBasis::Field SpectralBasis::evaluate(
  const std::vector<double> &c,
  const Vec2 y) const
{
  Workspace ws;
  return evaluate(c, y, ws);
}

}  // namespace conformality
