#include "conformality/mesh_check.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace conformality
{

std::vector<MeshAreaCheck> check_mesh_areas(
  const SpectralBasis &basis,
  const std::vector<AcceptedStep> &steps,
  const double lx,
  const double ly,
  const std::size_t base_cells,
  const std::size_t max_refinement,
  const ImplicitMidpointSettings &solver)
{
  if (
    base_cells == 0 || max_refinement == 0 ||
    (max_refinement & (max_refinement - 1)) != 0) {
    throw std::invalid_argument(
      "check_mesh_areas: base_cells must be positive and max_refinement a "
      "power of 2");
  }
  const double longer = std::max(lx, ly);
  const auto ncx = std::max<std::size_t>(
    1,
    static_cast<std::size_t>(
      std::lround(static_cast<double>(base_cells) * lx / longer)));
  const auto ncy = std::max<std::size_t>(
    1,
    static_cast<std::size_t>(
      std::lround(static_cast<double>(base_cells) * ly / longer)));
  const std::size_t nx = ncx * max_refinement;
  const std::size_t ny = ncy * max_refinement;
  const double fx = lx / static_cast<double>(nx);
  const double fy = ly / static_cast<double>(ny);
  std::vector<Vec2> pts;
  pts.reserve((nx + 1) * (ny + 1));
  for (std::size_t i = 0; i <= nx; ++i) {
    for (std::size_t j = 0; j <= ny; ++j) {
      pts.push_back(
        {static_cast<double>(i) * fx, static_cast<double>(j) * fy});
    }
  }
  static_cast<void>(replay(pts, nullptr, basis, steps, solver));
  auto at = [&](const std::size_t i, const std::size_t j) {
    return pts[i * (ny + 1) + j];
  };

  // Mesh of cells spanning `stride` lattice intervals, with boundary
  // samples every `step` intervals
  auto evaluate = [&](const std::size_t stride, const std::size_t step) {
    MeshAreaCheck m;
    m.cells_long = base_cells * (max_refinement / stride);
    m.samples_per_edge = stride / step;
    m.n_cells = (nx / stride) * (ny / stride);
    const double nominal =
      static_cast<double>(stride) * fx * static_cast<double>(stride) * fy;
    double sum = 0.0;
    for (std::size_t i0 = 0; i0 < nx; i0 += stride) {
      for (std::size_t j0 = 0; j0 < ny; j0 += stride) {

        // Shoelace formula relative to the first corner, counterclockwise
        const Vec2 origin = at(i0, j0);
        double twice_area = 0.0;
        Vec2 prev{0.0, 0.0};
        auto visit = [&](const std::size_t i, const std::size_t j) {
          const Vec2 q = at(i, j);
          const Vec2 p{q.x - origin.x, q.y - origin.y};
          twice_area += prev.x * p.y - p.x * prev.y;
          prev = p;
        };
        for (std::size_t k = step; k <= stride; k += step) {
          visit(i0 + k, j0);
        }
        for (std::size_t k = step; k <= stride; k += step) {
          visit(i0 + stride, j0 + k);
        }
        for (std::size_t k = step; k <= stride; k += step) {
          visit(i0 + stride - k, j0 + stride);
        }
        for (std::size_t k = step; k <= stride; k += step) {
          visit(i0, j0 + stride - k);
        }
        const double rel = std::abs(0.5 * twice_area / nominal - 1.0);
        if (rel > m.max_rel_error) {
          m.max_rel_error = rel;
          m.worst_x = static_cast<double>(i0) * fx;
          m.worst_y = static_cast<double>(j0) * fy;
        }
        sum += rel;
      }
    }
    m.mean_rel_error = sum / static_cast<double>(m.n_cells);
    return m;
  };
  std::vector<MeshAreaCheck> checks;
  for (std::size_t r = 1; r <= max_refinement; r *= 2) {
    checks.push_back(evaluate(max_refinement / r, max_refinement / r));
  }
  for (std::size_t r = 2; r <= max_refinement; r *= 2) {
    checks.push_back(evaluate(max_refinement, max_refinement / r));
  }
  return checks;
}

}  // namespace conformality
