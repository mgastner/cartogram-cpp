#ifndef CONFORMALITY_TEST_MAPS_HPP_
#define CONFORMALITY_TEST_MAPS_HPP_

// Manufactured maps and helpers shared by the conformality tests

#include "conformality/node_set.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <random>
#include <vector>

namespace conformality_test
{

using conformality::det;
using conformality::make_cell_centre_nodes;
using conformality::NodeSet;

inline constexpr double pi = std::numbers::pi;

// Handover test map on [0, 1]^2: T(x) = (x1 + eps sin(2 pi x1) / (2 pi), x2).
// Its energy for alpha = 0 is (1 - eps^2)^(-1/2) - 1.
inline NodeSet shear_map_nodes(double eps, std::size_t n)
{
  NodeSet nodes = make_cell_centre_nodes(1.0, 1.0, n, n);
  for (std::size_t j = 0; j < nodes.size(); ++j) {
    const double x1 = nodes.x[j].x;
    nodes.y[j] = {x1 + eps * std::sin(2 * pi * x1) / (2 * pi), nodes.x[j].y};
    nodes.jac[j] = {1.0 + eps * std::cos(2 * pi * x1), 0.0, 0.0, 1.0};
    nodes.g0[j] = det(nodes.jac[j]);
  }
  return nodes;
}

// Side-preserving map on [0, lx] x [0, ly] without separable structure or
// mirror symmetry, so that every modal derivative is generically nonzero:
//   T1 = x1 + a sin(pi x1/lx) sin(pi x2/ly) (1 + x1/lx + 2 x2/ly)
//   T2 = x2 + b sin(2 pi x1/lx) sin(pi x2/ly) (1 + x1/lx)
inline NodeSet asymmetric_map_nodes(
  double lx,
  double ly,
  std::size_t nx,
  std::size_t ny)
{
  constexpr double a = 0.05;
  constexpr double b = 0.03;
  NodeSet nodes = make_cell_centre_nodes(lx, ly, nx, ny);
  for (std::size_t j = 0; j < nodes.size(); ++j) {
    const double x1 = nodes.x[j].x;
    const double x2 = nodes.x[j].y;
    const double s1 = std::sin(pi * x1 / lx), c1 = std::cos(pi * x1 / lx);
    const double s2 = std::sin(pi * x2 / ly), c2 = std::cos(pi * x2 / ly);
    const double s1d = std::sin(2 * pi * x1 / lx);
    const double c1d = std::cos(2 * pi * x1 / lx);
    const double psi = 1.0 + x1 / lx + 2.0 * x2 / ly;
    const double chi = 1.0 + x1 / lx;
    nodes.y[j] = {x1 + a * s1 * s2 * psi, x2 + b * s1d * s2 * chi};
    nodes.jac[j] = {
      1.0 + a * ((pi / lx) * c1 * s2 * psi + s1 * s2 / lx),
      a * ((pi / ly) * s1 * c2 * psi + 2.0 * s1 * s2 / ly),
      b * ((2 * pi / lx) * c1d * s2 * chi + s1d * s2 / lx),
      1.0 + b * (pi / ly) * s1d * c2 * chi};
    nodes.g0[j] = det(nodes.jac[j]);
  }
  return nodes;
}

inline double max_abs(const std::vector<double> &v)
{
  double m = 0.0;
  for (const double x : v) {
    m = std::max(m, std::abs(x));
  }
  return m;
}

inline std::vector<double> random_coefficients(
  std::size_t n,
  double amplitude,
  std::mt19937 &rng)
{
  std::uniform_real_distribution<double> dist(-amplitude, amplitude);
  std::vector<double> c(n);
  for (auto &x : c) {
    x = dist(rng);
  }
  return c;
}

}  // namespace conformality_test

#endif  // CONFORMALITY_TEST_MAPS_HPP_
