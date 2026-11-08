// Replays a logged conformality correction A and checks, independently of
// the run that produced it, that A preserves area: a regular mesh of the
// final lattice is mapped by A, and the areas of the image cells are
// compared with the cells' areas (see conformality/mesh_check.hpp).
//
// Usage: conformality_check STEPS_JSON [--base N] [--refine R] [--threads T]
//
// STEPS_JSON is the <inset>_conformality_steps.json file that cartogram
// writes with --conf_optimize. The mesh has N cells along the longer side
// (default 64); the checks cover N, 2N, ..., N R cells with one sample per
// edge and N cells with 2, 4, ..., R samples per edge (R a power of 2,
// default 8). The replay uses N R intervals along the longer side.

#include "conformality/mesh_check.hpp"
#include "conformality/optimizer.hpp"
#include "conformality/parallel.hpp"
#include "conformality/spectral_basis.hpp"
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <utility>
#include <vector>

namespace
{

int usage()
{
  std::cerr << "usage: conformality_check STEPS_JSON [--base N] "
               "[--refine R] [--threads T]"
            << std::endl;
  return 2;
}

}  // namespace

int main(int argc, char *argv[])
{
  if (argc < 2) {
    return usage();
  }
  std::size_t base = 64;
  std::size_t refine = 8;
  int threads = 0;
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    if (i + 1 >= argc) {
      return usage();
    }
    const auto value = std::strtoull(argv[++i], nullptr, 10);
    if (arg == "--base") {
      base = static_cast<std::size_t>(value);
    } else if (arg == "--refine") {
      refine = static_cast<std::size_t>(value);
    } else if (arg == "--threads") {
      threads = static_cast<int>(value);
    } else {
      return usage();
    }
  }
  if (threads > 0) {
    conformality::set_max_threads(threads);
  }

  try {
    std::ifstream in(argv[1]);
    if (!in) {
      std::cerr << "cannot read " << argv[1] << std::endl;
      return 1;
    }
    const nlohmann::json log = nlohmann::json::parse(in);
    const double lx = log.at("lattice").at(0).get<double>();
    const double ly = log.at("lattice").at(1).get<double>();
    const auto mx = log.at("modes").at(0).get<std::size_t>();
    const auto my = log.at("modes").at(1).get<std::size_t>();
    conformality::FlowSettings solver;
    // Logs written before the integrator option used implicit midpoint
    solver.integrator = conformality::integrator_from_string(
      log.at("solver").value("integrator", std::string("midpoint")));
    solver.residual_tol = log.at("solver").at("residual_tol").get<double>();
    solver.max_iterations =
      log.at("solver").at("max_iterations").get<unsigned>();
    std::vector<conformality::AcceptedStep> steps;
    for (const auto &s : log.at("steps")) {
      conformality::AcceptedStep step;
      step.coeffs = s.at("coeffs").get<std::vector<double>>();
      step.h = s.at("h").get<double>();
      step.n_substeps = s.at("n_substeps").get<std::size_t>();
      steps.push_back(std::move(step));
    }
    const conformality::SpectralBasis basis(lx, ly, mx, my);
    std::cout << "lattice " << lx << " x " << ly << ", modes " << mx << " x "
              << my << ", " << steps.size() << " steps, "
              << conformality::max_threads() << " threads; replay of "
              << base * refine << " intervals along the longer side"
              << std::endl;

    const auto t0 = std::chrono::steady_clock::now();
    const auto checks = conformality::check_mesh_areas(
      basis,
      steps,
      lx,
      ly,
      base,
      refine,
      solver);
    const auto t1 = std::chrono::steady_clock::now();

    std::cout << std::setprecision(3) << std::scientific;
    std::cout << "cells_long samples_per_edge n_cells max_rel_error "
                 "ratio mean_rel_error worst_cell_corner"
              << std::endl;
    // Ratios compare each check with the previous check of its series. The
    // first check (base cells, one sample per edge) starts both series.
    double previous = 0.0;
    for (std::size_t k = 0; k < checks.size(); ++k) {
      const auto &c = checks[k];
      const bool starts_samples_series = c.samples_per_edge == 2;
      const double reference =
        starts_samples_series ? checks.front().max_rel_error : previous;
      std::cout << c.cells_long << ' ' << c.samples_per_edge << ' '
                << c.n_cells << ' ' << c.max_rel_error << ' ';
      if (k > 0 && reference > 0.0) {
        std::cout << std::fixed << std::setprecision(2)
                  << reference / c.max_rel_error << std::scientific
                  << std::setprecision(3);
      } else {
        std::cout << '-';
      }
      std::cout << ' ' << c.mean_rel_error << " (" << std::fixed
                << std::setprecision(1) << c.worst_x << ", " << c.worst_y
                << ")" << std::scientific << std::setprecision(3) << std::endl;
      previous = c.max_rel_error;
    }
    std::cout << std::fixed << std::setprecision(1) << "replay and check "
              << std::chrono::duration<double>(t1 - t0).count() << " s"
              << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
