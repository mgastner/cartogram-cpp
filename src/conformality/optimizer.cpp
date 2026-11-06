#include "conformality/optimizer.hpp"
#include "conformality/objective.hpp"
#include <algorithm>
#include <cmath>
#include <ostream>
#include <stdexcept>
#include <utility>

namespace conformality
{

namespace
{

void merge(FlowStats &total, const FlowStats &part)
{
  total.converged = total.converged && part.converged;
  total.max_iterations = std::max(total.max_iterations, part.max_iterations);
  total.max_residual = std::max(total.max_residual, part.max_residual);
}

}  // namespace

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
    !(settings.substep_factor > 0.0) || settings.max_substeps == 0) {
    throw std::invalid_argument("optimize: invalid settings");
  }

  OptimizationResult result;
  double energy = objective(nodes);
  result.initial_energy = energy;
  result.final_energy = energy;
  if (!std::isfinite(energy)) {
    result.outcome = Outcome::numerical_failure;
    result.reason = "initial map has a non-positive Jacobian determinant";
    return result;
  }

  double h = settings.initial_step;
  double r0 = 0.0;
  for (unsigned k = 0;; ++k) {
    const std::vector<double> f = modal_derivatives(nodes, basis);
    const H1Descent descent = h1_descent(basis, f);

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
    h = std::min(h, settings.max_step);
    const double bound = velocity_gradient_bound(basis, descent.coeffs);

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
                                    ? objective(trial)
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
    result.history.push_back(record);
    if (!accepted) {
      result.outcome = Outcome::numerical_failure;
      result.reason = "no trial satisfied the acceptance tests";
      break;
    }

    nodes = std::move(trial);
    const double previous = energy;
    energy = record.energy_after;
    result.final_energy = energy;
    if (previous - energy <= settings.energy_tol * previous) {
      result.outcome = Outcome::converged;
      result.reason = "relative energy decrease below tolerance";
      break;
    }
    h = std::min(2.0 * h, settings.max_step);
  }
  return result;
}

FlowStats replay(
  NodeSet &nodes,
  const SpectralBasis &basis,
  const std::vector<AcceptedStep> &steps,
  const ImplicitMidpointSettings &solver)
{
  FlowStats total;
  for (const auto &step : steps) {
    merge(
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
  const ImplicitMidpointSettings &solver)
{
  FlowStats total;
  for (const auto &step : steps) {
    merge(
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
