#ifndef CONFORMALITY_MESH_CHECK_HPP_
#define CONFORMALITY_MESH_CHECK_HPP_

#include "conformality/flow.hpp"
#include "conformality/optimizer.hpp"
#include "conformality/spectral_basis.hpp"
#include <cstddef>
#include <vector>

namespace conformality
{

// Position-based check of area preservation by a replayed correction A on
// [0, lx] x [0, ly]. A regular mesh of cells of size h is mapped by A; each
// image cell is the polygon through s samples per edge. Besides any loss of
// area preservation by A, the polygon misses the curvature of the images of
// the edges, an O((h / s)^2) error. Two refinements separate the two:
//   - finer meshes with s = 1: the error falls by about 4 per halving of h,
//     once h is below the length scale on which DA varies;
//   - more samples per edge on a fixed mesh: the error tends to the area
//     error of A itself.
struct MeshAreaCheck {
  std::size_t cells_long{0};
  std::size_t samples_per_edge{1};
  std::size_t n_cells{0};
  double max_rel_error{0.0};
  double mean_rel_error{0.0};

  // Lower-left corner, before A, of the cell with the largest error
  double worst_x{0.0};
  double worst_y{0.0};
};

// Replays A once on a lattice with base_cells * max_refinement intervals
// along the longer side (max_refinement a power of 2) and returns the
// checks for meshes of base_cells * r cells with s = 1, followed by the
// mesh of base_cells cells with s = r, for r = 1, 2, 4, ..., max_refinement
[[nodiscard]] std::vector<MeshAreaCheck> check_mesh_areas(
  const SpectralBasis &basis,
  const std::vector<AcceptedStep> &steps,
  double lx,
  double ly,
  std::size_t base_cells,
  std::size_t max_refinement,
  const ImplicitMidpointSettings &solver = {});

}  // namespace conformality

#endif  // CONFORMALITY_MESH_CHECK_HPP_
