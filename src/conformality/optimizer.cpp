#include "conformality/optimizer.hpp"
#include "conformality/objective.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ostream>
#include <random>
#include <stdexcept>
#include <utility>

namespace conformality
{

std::string to_string(const Outcome outcome)
{
  switch (outcome) {
  case Outcome::converged:
    return "converged";
  case Outcome::numerical_failure:
    return "numerical_failure";
  case Outcome::budget_exhausted:
    return "budget_exhausted";
  }
  return "unknown";
}

std::ostream &operator<<(std::ostream &os, const Outcome outcome)
{
  return os << to_string(outcome);
}

double max_det_drift(const NodeSet &nodes)
{
  double drift = 0.0;
  for (std::size_t j = 0; j < nodes.size(); ++j) {
    drift = std::max(drift, std::abs(det(nodes.jac[j]) / nodes.g0[j] - 1.0));
  }
  return drift;
}

OptimizationResult optimize(
  NodeSet &nodes,
  const SpectralBasis &basis,
  const OptimizerSettings &settings)
{
  nodes.check_consistent();
  if (
    !(settings.sigma > 0.0 && settings.sigma < 1.0) ||
    !(settings.substep_factor > 0.0) || settings.max_substeps == 0 ||
    !(settings.max_strain > 0.0) || !(settings.stall_tol >= 0.0) ||
    static_cast<bool>(settings.energy) !=
      static_cast<bool>(settings.derivatives)) {
    throw std::invalid_argument("optimize: invalid settings");
  }
  const auto energy_of = [&settings](const NodeSet &n) {
    return settings.energy ? settings.energy(n) : objective(n);
  };
  const auto derivatives_of = [&settings, &basis](const NodeSet &n) {
    return settings.derivatives ? settings.derivatives(n)
                                : modal_derivatives(n, basis);
  };

  OptimizationResult result;
  double energy = energy_of(nodes);
  result.initial_energy = energy;
  result.final_energy = energy;
  if (!std::isfinite(energy)) {
    result.outcome = Outcome::numerical_failure;
    result.reason = "initial map has a non-positive Jacobian determinant";
    return result;
  }

  // Inverse H1 metric, indexed like the coefficients
  std::vector<double> inv_g(basis.n_modes());
  for (std::size_t m = 1; m <= basis.mx(); ++m) {
    for (std::size_t n = 1; n <= basis.my(); ++n) {
      inv_g[basis.index(m, n)] = 1.0 / basis.h1_norm_sq(m, n);
    }
  }
  auto dot = [](const std::vector<double> &a, const std::vector<double> &b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
      sum += a[i] * b[i];
    }
    return sum;
  };

  // L-BFGS memory: steps s_i, gradient changes y_i and 1 / (s_i . y_i)
  std::vector<std::vector<double>> mem_s;
  std::vector<std::vector<double>> mem_y;
  std::vector<double> mem_rho;
  std::vector<double> previous_f;
  std::vector<double> previous_s;

  double h = settings.initial_step;
  double r0 = 0.0;
  for (unsigned k = 0;; ++k) {
    const auto iteration_start = std::chrono::steady_clock::now();
    const std::vector<double> f = derivatives_of(nodes);
    H1Descent descent = h1_descent(basis, f);

    // Quasi-Newton direction from the stored pairs
    bool quasi_newton = false;
    if (settings.lbfgs_memory > 0) {
      if (!previous_s.empty()) {
        std::vector<double> y(f.size());
        for (std::size_t i = 0; i < f.size(); ++i) {
          y[i] = f[i] - previous_f[i];
        }
        const double sy = dot(previous_s, y);
        if (sy > 0.0) {
          if (mem_s.size() == settings.lbfgs_memory) {
            mem_s.erase(mem_s.begin());
            mem_y.erase(mem_y.begin());
            mem_rho.erase(mem_rho.begin());
          }
          mem_s.push_back(previous_s);
          mem_y.push_back(std::move(y));
          mem_rho.push_back(1.0 / sy);
        }
      }
      if (!mem_s.empty()) {
        std::vector<double> q = f;
        std::vector<double> a(mem_s.size());
        for (std::size_t i = mem_s.size(); i-- > 0;) {
          a[i] = mem_rho[i] * dot(mem_s[i], q);
          for (std::size_t j = 0; j < q.size(); ++j) {
            q[j] -= a[i] * mem_y[i][j];
          }
        }
        const std::vector<double> &yn = mem_y.back();
        double y_ginv_y = 0.0;
        for (std::size_t j = 0; j < yn.size(); ++j) {
          y_ginv_y += yn[j] * yn[j] * inv_g[j];
        }
        const double gamma = 1.0 / (mem_rho.back() * y_ginv_y);
        for (std::size_t j = 0; j < q.size(); ++j) {
          q[j] *= gamma * inv_g[j];
        }
        for (std::size_t i = 0; i < mem_s.size(); ++i) {
          const double b = mem_rho[i] * dot(mem_y[i], q);
          for (std::size_t j = 0; j < q.size(); ++j) {
            q[j] += (a[i] - b) * mem_s[i][j];
          }
        }
        std::vector<double> p(q.size());
        for (std::size_t j = 0; j < q.size(); ++j) {
          p[j] = -q[j];
        }
        if (dot(p, f) < 0.0) {
          descent.coeffs = std::move(p);
          quasi_newton = true;
        } else {
          mem_s.clear();
          mem_y.clear();
          mem_rho.clear();
        }
      }
    }

    // Slope of E along the velocity that will actually be integrated
    double slope = 0.0;
    for (std::size_t i = 0; i < f.size(); ++i) {
      slope += descent.coeffs[i] * f[i];
    }
    if (k == 0) {
      r0 = descent.r;
      result.initial_r = r0;
    }
    result.final_r = descent.r;

    if (!(r0 > 0.0) || std::sqrt(descent.r / r0) <= settings.gradient_tol) {
      result.outcome = Outcome::converged;
      result.reason = "projected gradient norm below tolerance";
      break;
    }
    if (k >= settings.max_iterations) {
      result.outcome = Outcome::budget_exhausted;
      result.reason = "maximum number of iterations reached";
      break;
    }
    if (!(slope < 0.0)) {
      result.outcome = Outcome::numerical_failure;
      result.reason = "velocity is not a descent direction";
      break;
    }

    IterationRecord record;
    record.iteration = k;
    record.energy = energy;
    record.r = descent.r;
    record.slope = slope;

    if (!(h > 0.0)) {
      h = energy / descent.r;
    }
    // Steepest-descent duration, kept for iterations without a quasi-Newton
    // direction
    const double h_descent = h;
    if (quasi_newton) {
      h = 1.0;
    }
    h = std::min(h, settings.max_step);
    const double bound = velocity_gradient_bound(basis, descent.coeffs);
    if (bound > 0.0) {
      h = std::min(h, settings.max_strain / bound);
    }

    bool accepted = false;
    NodeSet trial;
    for (unsigned b = 0; b <= settings.max_backtracks; ++b) {
      double n_real = std::ceil(h * bound / settings.substep_factor);
      if (n_real > static_cast<double>(settings.max_substeps)) {
        h = static_cast<double>(settings.max_substeps) *
            settings.substep_factor / bound;
        n_real = static_cast<double>(settings.max_substeps);
      }
      const std::size_t n_substeps =
        std::max<std::size_t>(1, static_cast<std::size_t>(n_real));

      trial = nodes;
      const FlowStats flow = advance_nodes(
        trial,
        basis,
        descent.coeffs,
        h,
        n_substeps,
        settings.solver);
      const double trial_energy = flow.converged
                                    ? energy_of(trial)
                                    : std::numeric_limits<double>::infinity();
      const double drift = max_det_drift(trial);
      const double margin = energy + settings.sigma * h * slope - trial_energy;

      record.backtracks = b;
      record.max_solver_iterations = flow.max_iterations;
      record.max_solver_residual = flow.max_residual;
      if (
        std::isfinite(trial_energy) && margin >= 0.0 &&
        drift <= settings.det_drift_tol) {
        record.h = h;
        record.strain = h * bound;
        record.n_substeps = n_substeps;
        record.energy_after = trial_energy;
        record.armijo_margin = margin;
        record.max_det_drift = drift;
        result.steps.push_back({descent.coeffs, h, n_substeps});
        accepted = true;
        break;
      }
      h *= 0.5;
    }
    record.seconds = std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - iteration_start)
                       .count();
    result.history.push_back(record);
    if (settings.on_iteration) {
      settings.on_iteration(record);
    }
    if (!accepted) {
      result.outcome = Outcome::numerical_failure;
      result.reason = "no trial satisfied the acceptance tests";
      break;
    }

    nodes = std::move(trial);
    if (settings.lbfgs_memory > 0) {
      previous_f = f;
      previous_s = descent.coeffs;
      for (double &x : previous_s) {
        x *= h;
      }
    }
    const double previous = energy;
    energy = record.energy_after;
    result.final_energy = energy;
    if (previous - energy <= settings.energy_tol * previous) {
      result.outcome = Outcome::converged;
      result.reason = "relative energy decrease below tolerance";
      break;
    }
    const std::size_t window = settings.stall_window;
    if (window > 0 && result.history.size() >= window) {
      const double earlier =
        result.history[result.history.size() - window].energy;
      if (earlier - energy <= settings.stall_tol * energy) {
        result.outcome = Outcome::converged;
        result.reason = "relative energy decrease over the last " +
                        std::to_string(window) + " iterations below tolerance";
        break;
      }
    }
    // The next steepest-descent trial doubles the accepted duration; after a
    // quasi-Newton step, it reuses the steepest-descent duration of this
    // iteration
    h = quasi_newton ? h_descent : std::min(2.0 * h, settings.max_step);
  }
  return result;
}

AcceptedStep random_step(
  const SpectralBasis &basis,
  const double strain,
  const std::uint64_t seed,
  const double substep_factor)
{
  if (
    !(strain > 0.0) || !std::isfinite(strain) || !(substep_factor > 0.0) ||
    !std::isfinite(substep_factor)) {
    throw std::invalid_argument(
      "random_step: strain and substep_factor must be positive and finite");
  }
  std::mt19937_64 rng(seed);
  // Uniform deviate in (0, 1] from the top 53 bits of a draw
  auto uniform = [&rng]() {
    return (static_cast<double>(rng() >> 11) + 1.0) * 0x1.0p-53;
  };
  constexpr double two_pi = 6.283185307179586;
  AcceptedStep step;
  step.coeffs.assign(basis.n_modes(), 0.0);
  for (std::size_t m = 1; m <= basis.mx(); ++m) {
    for (std::size_t n = 1; n <= basis.my(); ++n) {
      const double u1 = uniform();
      const double u2 = uniform();
      const double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(two_pi * u2);
      const double k2 = basis.kx(m) * basis.kx(m) + basis.ky(n) * basis.ky(n);
      step.coeffs[basis.index(m, n)] = z / (k2 * std::sqrt(k2));
    }
  }
  step.h = strain / velocity_gradient_bound(basis, step.coeffs);
  step.n_substeps =
    static_cast<std::size_t>(std::ceil(strain / substep_factor));
  return step;
}

FlowStats replay(
  NodeSet &nodes,
  const SpectralBasis &basis,
  const std::vector<AcceptedStep> &steps,
  const FlowSettings &solver)
{
  FlowStats total;
  for (const auto &step : steps) {
    accumulate(
      total,
      advance_nodes(
        nodes,
        basis,
        step.coeffs,
        step.h,
        step.n_substeps,
        solver));
  }
  return total;
}

FlowStats replay(
  std::vector<Vec2> &points,
  std::vector<Mat2> *tangents,
  const SpectralBasis &basis,
  const std::vector<AcceptedStep> &steps,
  const FlowSettings &solver)
{
  FlowStats total;
  for (const auto &step : steps) {
    accumulate(
      total,
      advance_points(
        points,
        tangents,
        basis,
        step.coeffs,
        step.h,
        step.n_substeps,
        solver));
  }
  return total;
}

}  // namespace conformality
