#ifndef CONFORMALITY_SPECTRAL_BASIS_HPP_
#define CONFORMALITY_SPECTRAL_BASIS_HPP_

#include "conformality/types.hpp"
#include <cstddef>
#include <vector>

namespace conformality
{

// Fills sin_out[j - 1] = sin(j * theta) and cos_out[j - 1] = cos(j * theta)
// for j = 1, ..., count. Uses the angle-addition recurrence, whose rounding
// error grows only linearly in j.
void fill_harmonics(
  double theta,
  std::size_t count,
  double *sin_out,
  double *cos_out) noexcept;

// Fourier sine basis for stream functions on the rectangle
// [0, lx] x [0, ly]:
//
//   s(y) = sum_{m=1}^{mx} sum_{n=1}^{my} c_mn sin(kx_m y1) sin(ky_n y2),
//   kx_m = m pi / lx,  ky_n = n pi / ly.
//
// The velocity v = (ds/dy2, -ds/dy1) is divergence free and tangent to every
// side of the rectangle, and it vanishes at the corners. Coefficients are
// stored in a vector of length mx * my, with c_mn at index(m, n).
class SpectralBasis
{
public:
  // Value of the stream function, the velocity, and the velocity gradient
  // (dv)_ij = dv_i/dy_j at one point
  struct Field {
    double s{0.0};
    Vec2 v{};
    Mat2 dv{};
  };

  // Scratch space for evaluate(). Reusing one workspace across calls avoids
  // reallocating the harmonic tables.
  struct Workspace {
    std::vector<double> sin_x, cos_x, sin_y, cos_y;
  };

  // Throws std::invalid_argument unless lx > 0, ly > 0, mx >= 1, my >= 1
  SpectralBasis(double lx, double ly, std::size_t mx, std::size_t my);

  [[nodiscard]] double lx() const noexcept;
  [[nodiscard]] double ly() const noexcept;
  [[nodiscard]] std::size_t mx() const noexcept;
  [[nodiscard]] std::size_t my() const noexcept;
  [[nodiscard]] std::size_t n_modes() const noexcept;

  // Position of c_mn in the coefficient vector (1 <= m <= mx, 1 <= n <= my)
  [[nodiscard]] std::size_t index(std::size_t m, std::size_t n) const noexcept;

  [[nodiscard]] double kx(std::size_t m) const noexcept;
  [[nodiscard]] double ky(std::size_t n) const noexcept;

  // Squared homogeneous H1 norm of the velocity mode (m, n),
  // G_mn = integral of dv_mn : dv_mn = (lx ly / 4) (kx_m^2 + ky_n^2)^2
  [[nodiscard]] double h1_norm_sq(std::size_t m, std::size_t n) const noexcept;

  // Evaluates s, v, and dv at y by direct summation over all modes.
  // Precondition: c.size() == n_modes().
  [[nodiscard]] Field evaluate(
    const std::vector<double> &c,
    Vec2 y,
    Workspace &ws) const;

  // Convenience overload that allocates its own workspace
  [[nodiscard]] Field evaluate(const std::vector<double> &c, Vec2 y) const;

private:
  double lx_;
  double ly_;
  std::size_t mx_;
  std::size_t my_;

  // kx_[m - 1] = kx_m, kx_sq_[m - 1] = kx_m^2, and likewise for y
  std::vector<double> kx_, kx_sq_, ky_, ky_sq_;
};

}  // namespace conformality

#endif  // CONFORMALITY_SPECTRAL_BASIS_HPP_
