#include "conformality/flow.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace conformality
{

namespace
{

double distance(const Vec2 a, const Vec2 b) noexcept
{
  return std::hypot(a.x - b.x, a.y - b.y);
}

void check_flow_arguments(const double h, const std::size_t n_substeps)
{
  if (!std::isfinite(h) || n_substeps == 0) {
    throw std::invalid_argument(
      "advance: h must be finite and n_substeps positive");
  }
}

// Advances one point (and optionally a tangent) through n_substeps
void advance_one(
  Vec2 &y,
  Mat2 *tangent,
  const SpectralBasis &basis,
  const std::vector<double> &c,
  const double dt,
  const std::size_t n_substeps,
  const ImplicitMidpointSettings &settings,
  SpectralBasis::Workspace &ws,
  FlowStats &stats)
{
  for (std::size_t step = 0; step < n_substeps; ++step) {
    const SubstepResult r =
      implicit_midpoint_substep(basis, c, y, dt, settings, ws);
    y = r.y;
    if (tangent != nullptr) {
      *tangent = r.tangent * (*tangent);
    }
    stats.converged = stats.converged && r.converged;
    stats.max_iterations = std::max(stats.max_iterations, r.iterations);
    stats.max_residual = std::max(stats.max_residual, r.residual);
  }
}

// Advances y[j] (and tangent[j] if tangent is not null) for all
// j < count. Every point is integrated on its own, so neither the thread
// count nor the schedule affects the result.
FlowStats advance_all(
  Vec2 *y,
  Mat2 *tangent,
  const std::size_t count,
  const SpectralBasis &basis,
  const std::vector<double> &c,
  const double dt,
  const std::size_t n_substeps,
  const ImplicitMidpointSettings &settings)
{
  FlowStats stats;
#pragma omp parallel
  {
    SpectralBasis::Workspace ws;
    FlowStats local;

    // Fixed-point iteration counts vary between points, so hand out small
    // chunks dynamically
#pragma omp for schedule(dynamic, 64) nowait
    for (std::size_t j = 0; j < count; ++j) {
      advance_one(
        y[j],
        tangent != nullptr ? tangent + j : nullptr,
        basis,
        c,
        dt,
        n_substeps,
        settings,
        ws,
        local);
    }
#pragma omp critical(conformality_flow_stats)
    accumulate(stats, local);
  }
  return stats;
}

}  // namespace

void accumulate(FlowStats &total, const FlowStats &part) noexcept
{
  total.converged = total.converged && part.converged;
  total.max_iterations = std::max(total.max_iterations, part.max_iterations);
  total.max_residual = std::max(total.max_residual, part.max_residual);
}

Mat2 cayley(const Mat2 &d, const double dt) noexcept
{
  const double s = 0.5 * dt;

  // m = I - s D and its inverse via the adjugate
  const Mat2 m{1.0 - s * d.a11, -s * d.a12, -s * d.a21, 1.0 - s * d.a22};
  const double det_m = det(m);
  const Mat2 m_inv{
    m.a22 / det_m,
    -m.a12 / det_m,
    -m.a21 / det_m,
    m.a11 / det_m};
  const Mat2 p{1.0 + s * d.a11, s * d.a12, s * d.a21, 1.0 + s * d.a22};
  return m_inv * p;
}

double velocity_gradient_bound(
  const SpectralBasis &basis,
  const std::vector<double> &c)
{
  assert(c.size() == basis.n_modes());
  double bound = 0.0;
  for (std::size_t m = 1; m <= basis.mx(); ++m) {
    for (std::size_t n = 1; n <= basis.my(); ++n) {
      const double kx = basis.kx(m);
      const double ky = basis.ky(n);
      bound += std::abs(c[basis.index(m, n)]) * (kx * kx + ky * ky);
    }
  }
  return bound;
}

SubstepResult implicit_midpoint_substep(
  const SpectralBasis &basis,
  const std::vector<double> &c,
  const Vec2 y,
  const double dt,
  const ImplicitMidpointSettings &settings,
  SpectralBasis::Workspace &ws)
{
  const double tol = settings.residual_tol * std::max(basis.lx(), basis.ly());

  // Explicit Euler predictor, then fixed-point iteration on
  // y+ = y + dt v((y + y+) / 2)
  SubstepResult result;
  const Vec2 v0 = basis.evaluate(c, y, ws).v;
  Vec2 yp{y.x + dt * v0.x, y.y + dt * v0.y};
  for (unsigned it = 1; it <= settings.max_iterations; ++it) {
    const Vec2 mid{0.5 * (y.x + yp.x), 0.5 * (y.y + yp.y)};
    const Vec2 v = basis.evaluate(c, mid, ws).v;
    const Vec2 next{y.x + dt * v.x, y.y + dt * v.y};
    const double change = distance(next, yp);
    yp = next;
    result.iterations = it;
    if (!std::isfinite(change)) {
      break;
    }
    if (change <= tol) {
      result.converged = true;
      break;
    }
  }

  // Tangent and residual at the midpoint of the returned step
  const Vec2 mid{0.5 * (y.x + yp.x), 0.5 * (y.y + yp.y)};
  const auto f = basis.evaluate(c, mid, ws);
  result.y = yp;
  result.tangent = cayley(f.dv, dt);
  result.residual = distance({y.x + dt * f.v.x, y.y + dt * f.v.y}, yp);
  if (!std::isfinite(result.residual)) {
    result.converged = false;
  }
  return result;
}

FlowStats advance_nodes(
  NodeSet &nodes,
  const SpectralBasis &basis,
  const std::vector<double> &c,
  const double h,
  const std::size_t n_substeps,
  const ImplicitMidpointSettings &settings)
{
  nodes.check_consistent();
  check_flow_arguments(h, n_substeps);
  const double dt = h / static_cast<double>(n_substeps);
  return advance_all(
    nodes.y.data(),
    nodes.jac.data(),
    nodes.size(),
    basis,
    c,
    dt,
    n_substeps,
    settings);
}

FlowStats advance_points(
  std::vector<Vec2> &points,
  std::vector<Mat2> *tangents,
  const SpectralBasis &basis,
  const std::vector<double> &c,
  const double h,
  const std::size_t n_substeps,
  const ImplicitMidpointSettings &settings)
{
  check_flow_arguments(h, n_substeps);
  if (tangents != nullptr && tangents->size() != points.size()) {
    throw std::invalid_argument(
      "advance_points: tangents and points differ in length");
  }
  const double dt = h / static_cast<double>(n_substeps);
  return advance_all(
    points.data(),
    tangents != nullptr ? tangents->data() : nullptr,
    points.size(),
    basis,
    c,
    dt,
    n_substeps,
    settings);
}

}  // namespace conformality
