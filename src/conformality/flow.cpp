#include "conformality/flow.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numbers>
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

// A Hamiltonian H(q, p) = sum_ij h[i * np + j] sin(kq_i q) sin(kp_j p) on
// [0, lq] x [0, lp]. The stream function s(y1, y2) is H with q = y1, p = y2
// and h = c. The adjoint symplectic Euler step of s is the ordinary step of
// H(q, p) = -s(p, q) with q = y2, p = y1, whose coefficients are -c^T.
struct SeparableHamiltonian {
  const double *h;
  std::size_t nq;
  std::size_t np;
  const double *kq;
  const double *kq_sq;
  const double *kp;
  const double *kp_sq;
  double lq;
  double lp;
};

SeparableHamiltonian stream_function(
  const SpectralBasis &basis,
  const std::vector<double> &c) noexcept
{
  return {
    c.data(),
    basis.mx(),
    basis.my(),
    basis.kx_all().data(),
    basis.kx_sq_all().data(),
    basis.ky_all().data(),
    basis.ky_sq_all().data(),
    basis.lx(),
    basis.ly()};
}

// -c^T, laid out as the coefficients of H(q, p) = -s(p, q)
void negated_transpose(
  const SpectralBasis &basis,
  const std::vector<double> &c,
  std::vector<double> &out)
{
  const std::size_t mx = basis.mx();
  const std::size_t my = basis.my();
  out.resize(mx * my);
  for (std::size_t m = 0; m < mx; ++m) {
    for (std::size_t n = 0; n < my; ++n) {
      out[n * mx + m] = -c[m * my + n];
    }
  }
}

SeparableHamiltonian swapped_stream_function(
  const SpectralBasis &basis,
  const std::vector<double> &neg_ct) noexcept
{
  return {
    neg_ct.data(),
    basis.my(),
    basis.mx(),
    basis.ky_all().data(),
    basis.ky_sq_all().data(),
    basis.kx_all().data(),
    basis.kx_sq_all().data(),
    basis.ly(),
    basis.lx()};
}

// One symplectic Euler step of H in (q, p), implicit in p:
//   p+ = p - dt H_q(q, p+),  q+ = q + dt H_p(q, p+).
// Returns y = (q+, p+) and the tangent with rows (q+, p+) and columns (q, p).
//
// With q fixed, H and its derivatives are sums over j with the coefficients
//   a_j = sum_i h_ij kq_i cos(kq_i q),  b_j = sum_i h_ij sin(kq_i q),
//   e_j = sum_i h_ij kq_i^2 sin(kq_i q),
// and with T_j = sin(kp_j p), U_j = cos(kp_j p):
//   H_q = sum a_j T_j,             H_qp = sum a_j kp_j U_j,
//   H_p = sum b_j kp_j U_j,        H_pp = -sum b_j kp_j^2 T_j,
//   H_qq = -sum e_j T_j.
// The scalar equation f(z) = z - p + dt H_q(q, z) = 0 is solved by Newton's
// method from the explicit predictor z = p - dt H_q(q, p). Differentiating the
// two equations gives, with alpha = 1 + dt H_qp at (q, p+),
//   dp+/dq = -dt H_qq / alpha,  dp+/dp = 1 / alpha,
//   dq+/dq = alpha - dt^2 H_pp H_qq / alpha,  dq+/dp = dt H_pp / alpha,
// whose determinant is 1. For fixed q, z -> z + dt H_q(q, z) is increasing
// when alpha > 0, and for fixed p+, q -> q+ is increasing when alpha > 0, so
// the step is one to one where alpha > 0; elsewhere it is reported as failed.
SubstepResult euler_step(
  const SeparableHamiltonian &ham,
  const double q,
  const double p,
  const double dt,
  const double tol,
  const unsigned max_iterations,
  SpectralBasis::Workspace &ws)
{
  ws.sin_x.resize(ham.nq);
  ws.cos_x.resize(ham.nq);
  ws.sin_y.resize(ham.np);
  ws.cos_y.resize(ham.np);
  ws.line_a.assign(ham.np, 0.0);
  ws.line_b.assign(ham.np, 0.0);
  ws.line_e.assign(ham.np, 0.0);
  double *sin_q = ws.sin_x.data();
  double *cos_q = ws.cos_x.data();
  double *sin_p = ws.sin_y.data();
  double *cos_p = ws.cos_y.data();
  double *a = ws.line_a.data();
  double *b = ws.line_b.data();
  double *e = ws.line_e.data();

  fill_harmonics(std::numbers::pi * q / ham.lq, ham.nq, sin_q, cos_q);
  for (std::size_t i = 0; i < ham.nq; ++i) {
    const double si = sin_q[i];
    const double kci = ham.kq[i] * cos_q[i];
    const double k2si = ham.kq_sq[i] * sin_q[i];
    const double *row = ham.h + i * ham.np;
    for (std::size_t j = 0; j < ham.np; ++j) {
      a[j] += row[j] * kci;
      b[j] += row[j] * si;
      e[j] += row[j] * k2si;
    }
  }

  // H_q and H_qp at (q, z)
  auto line = [&](const double z, double &h_q, double &h_qp) {
    fill_harmonics(std::numbers::pi * z / ham.lp, ham.np, sin_p, cos_p);
    h_q = 0.0;
    h_qp = 0.0;
    for (std::size_t j = 0; j < ham.np; ++j) {
      h_q += a[j] * sin_p[j];
      h_qp += a[j] * ham.kp[j] * cos_p[j];
    }
  };

  SubstepResult result;
  double h_q = 0.0;
  double h_qp = 0.0;
  line(p, h_q, h_qp);
  double z = p - dt * h_q;
  for (unsigned it = 1; it <= max_iterations; ++it) {
    line(z, h_q, h_qp);
    const double change = (z - p + dt * h_q) / (1.0 + dt * h_qp);
    z -= change;
    result.iterations = it;
    if (!std::isfinite(change)) {
      break;
    }
    if (std::abs(change) <= tol) {
      result.converged = true;
      break;
    }
  }

  // Derivatives at (q, p+)
  fill_harmonics(std::numbers::pi * z / ham.lp, ham.np, sin_p, cos_p);
  double h_p = 0.0;
  double h_pp = 0.0;
  double h_qq = 0.0;
  h_q = 0.0;
  h_qp = 0.0;
  for (std::size_t j = 0; j < ham.np; ++j) {
    const double kp = ham.kp[j];
    h_q += a[j] * sin_p[j];
    h_qp += a[j] * kp * cos_p[j];
    h_p += b[j] * kp * cos_p[j];
    h_pp -= b[j] * ham.kp_sq[j] * sin_p[j];
    h_qq -= e[j] * sin_p[j];
  }
  const double alpha = 1.0 + dt * h_qp;
  result.y = {q + dt * h_p, z};
  result.tangent = {
    alpha - dt * dt * h_pp * h_qq / alpha,
    dt * h_pp / alpha,
    -dt * h_qq / alpha,
    1.0 / alpha};
  result.residual = std::abs(z - p + dt * h_q);
  const Mat2 &t = result.tangent;
  if (
    !(alpha > 0.0) || !std::isfinite(result.residual) ||
    !std::isfinite(result.y.x) || !std::isfinite(t.a11) ||
    !std::isfinite(t.a12) || !std::isfinite(t.a21) || !std::isfinite(t.a22)) {
    result.converged = false;
  }
  return result;
}

// Symplectic Euler implicit in y2 (q = y1, p = y2, H = s)
SubstepResult euler_implicit_in_y2(
  const SeparableHamiltonian &s,
  const Vec2 y,
  const double dt,
  const double tol,
  const unsigned max_iterations,
  SpectralBasis::Workspace &ws)
{
  return euler_step(s, y.x, y.y, dt, tol, max_iterations, ws);
}

// Symplectic Euler implicit in y1 (q = y2, p = y1, H = -s(p, q)); the result
// is mapped back to the order (y1, y2)
SubstepResult euler_implicit_in_y1(
  const SeparableHamiltonian &swapped,
  const Vec2 y,
  const double dt,
  const double tol,
  const unsigned max_iterations,
  SpectralBasis::Workspace &ws)
{
  SubstepResult r = euler_step(swapped, y.y, y.x, dt, tol, max_iterations, ws);
  r.y = {r.y.y, r.y.x};
  const Mat2 t = r.tangent;
  r.tangent = {t.a22, t.a21, t.a12, t.a11};
  return r;
}

SubstepResult verlet_step(
  const SeparableHamiltonian &s,
  const SeparableHamiltonian &swapped,
  const Vec2 y,
  const double dt,
  const double tol,
  const unsigned max_iterations,
  SpectralBasis::Workspace &ws)
{
  const SubstepResult first =
    euler_implicit_in_y2(s, y, 0.5 * dt, tol, max_iterations, ws);
  SubstepResult second =
    euler_implicit_in_y1(swapped, first.y, 0.5 * dt, tol, max_iterations, ws);
  second.tangent = second.tangent * first.tangent;
  second.iterations += first.iterations;
  second.residual = std::max(first.residual, second.residual);
  second.converged = first.converged && second.converged;
  return second;
}

// Everything a substep needs that is shared by all points of one advance
struct Stepper {
  const SpectralBasis &basis;
  const std::vector<double> &c;
  const FlowSettings &settings;
  SeparableHamiltonian s;
  SeparableHamiltonian swapped;
  double tol;

  SubstepResult step(
    const Vec2 y,
    const double dt,
    SpectralBasis::Workspace &ws) const
  {
    switch (settings.integrator) {
    case Integrator::symplectic_euler:
      return euler_implicit_in_y2(s, y, dt, tol, settings.max_iterations, ws);
    case Integrator::stormer_verlet:
      return verlet_step(s, swapped, y, dt, tol, settings.max_iterations, ws);
    case Integrator::implicit_midpoint:
      break;
    }
    return implicit_midpoint_substep(basis, c, y, dt, settings, ws);
  }
};

// Advances one point (and optionally a tangent) through n_substeps
void advance_one(
  Vec2 &y,
  Mat2 *tangent,
  const Stepper &stepper,
  const double dt,
  const std::size_t n_substeps,
  SpectralBasis::Workspace &ws,
  FlowStats &stats)
{
  for (std::size_t step = 0; step < n_substeps; ++step) {
    const SubstepResult r = stepper.step(y, dt, ws);
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
  const FlowSettings &settings)
{
  if (c.size() != basis.n_modes()) {
    throw std::invalid_argument("advance: wrong number of coefficients");
  }
  std::vector<double> neg_ct;
  if (settings.integrator == Integrator::stormer_verlet) {
    negated_transpose(basis, c, neg_ct);
  }
  const Stepper stepper{
    basis,
    c,
    settings,
    stream_function(basis, c),
    swapped_stream_function(basis, neg_ct),
    settings.residual_tol * std::max(basis.lx(), basis.ly())};

  FlowStats stats;
#pragma omp parallel
  {
    SpectralBasis::Workspace ws;
    FlowStats local;

    // Iteration counts vary between points, so hand out small chunks
    // dynamically
#pragma omp for schedule(dynamic, 64) nowait
    for (std::size_t j = 0; j < count; ++j) {
      advance_one(
        y[j],
        tangent != nullptr ? tangent + j : nullptr,
        stepper,
        dt,
        n_substeps,
        ws,
        local);
    }
#pragma omp critical(conformality_flow_stats)
    accumulate(stats, local);
  }
  return stats;
}

}  // namespace

std::string to_string(const Integrator integrator)
{
  switch (integrator) {
  case Integrator::implicit_midpoint:
    return "midpoint";
  case Integrator::symplectic_euler:
    return "euler";
  case Integrator::stormer_verlet:
    return "verlet";
  }
  return "unknown";
}

Integrator integrator_from_string(const std::string &name)
{
  for (const Integrator i :
       {Integrator::implicit_midpoint,
        Integrator::symplectic_euler,
        Integrator::stormer_verlet}) {
    if (name == to_string(i)) {
      return i;
    }
  }
  throw std::invalid_argument(
    "unknown integrator \"" + name + "\" (midpoint, euler or verlet)");
}

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
  const FlowSettings &settings,
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

SubstepResult symplectic_euler_substep(
  const SpectralBasis &basis,
  const std::vector<double> &c,
  const Vec2 y,
  const double dt,
  const FlowSettings &settings,
  SpectralBasis::Workspace &ws,
  const EulerVariant variant)
{
  assert(c.size() == basis.n_modes());
  const double tol = settings.residual_tol * std::max(basis.lx(), basis.ly());
  if (variant == EulerVariant::implicit_in_y2) {
    return euler_implicit_in_y2(
      stream_function(basis, c),
      y,
      dt,
      tol,
      settings.max_iterations,
      ws);
  }
  negated_transpose(basis, c, ws.transposed);
  return euler_implicit_in_y1(
    swapped_stream_function(basis, ws.transposed),
    y,
    dt,
    tol,
    settings.max_iterations,
    ws);
}

SubstepResult stormer_verlet_substep(
  const SpectralBasis &basis,
  const std::vector<double> &c,
  const Vec2 y,
  const double dt,
  const FlowSettings &settings,
  SpectralBasis::Workspace &ws)
{
  assert(c.size() == basis.n_modes());
  negated_transpose(basis, c, ws.transposed);
  return verlet_step(
    stream_function(basis, c),
    swapped_stream_function(basis, ws.transposed),
    y,
    dt,
    settings.residual_tol * std::max(basis.lx(), basis.ly()),
    settings.max_iterations,
    ws);
}

FlowStats advance_nodes(
  NodeSet &nodes,
  const SpectralBasis &basis,
  const std::vector<double> &c,
  const double h,
  const std::size_t n_substeps,
  const FlowSettings &settings)
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
  const FlowSettings &settings)
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
