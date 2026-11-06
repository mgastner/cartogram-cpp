#ifndef CONFORMALITY_TRACERS_HPP_
#define CONFORMALITY_TRACERS_HPP_

#include "conformality/node_set.hpp"
#include "conformality/types.hpp"
#include <cstddef>
#include <vector>

// Tracers that record the cumulative cartogram map T0 for the conformality
// postprocessor. They are seeded at fixed points x of the initial lattice
// [0, lx0] x [0, ly0] and moved through every successful flattening pass by
// the same piecewise-affine map that projects the polygons. Positions and
// Jacobians stay in the current lattice frame, so grid doublings
// (adjust_grid) scale them like the polygons.
struct ConformalityTracers {

  // Whether tracers were seeded for this inset
  bool active{false};

  // Initial lattice size and quadrature parameters
  unsigned int lx0{0};
  unsigned int ly0{0};
  unsigned int quad{0};
  unsigned int eval_quad{0};

  // Number of successful passes and of grid doublings since seeding
  unsigned int n_passes{0};
  unsigned int n_doublings{0};

  // Whether polygons were simplified after seeding. If so, the output
  // polygons are no longer exact images of the input polygons under T0.
  bool simplified_after_seeding{false};

  // Quadrature nodes used by the optimizer, and the index of the GeoDiv
  // that contains each node's original position (-1 if none)
  conformality::NodeSet nodes;
  std::vector<int> region;

  // Independent, finer evaluation cloud (never used by the optimizer)
  conformality::NodeSet eval_nodes;
  std::vector<int> eval_region;

  // Tracers on the sides of the rectangle (side 0: x = 0, 1: x = lx,
  // 2: y = 0, 3: y = ly) and at its corners, in the order (0, 0), (lx, 0),
  // (0, ly), (lx, ly)
  std::vector<conformality::Vec2> boundary;
  std::vector<unsigned char> boundary_side;
  std::vector<conformality::Vec2> corners;

  // Finite-difference companions: for node fd_parent[k], tracers seeded at
  // x + fd_delta e1 and x + fd_delta e2
  double fd_delta{1e-4};
  std::vector<std::size_t> fd_parent;
  std::vector<conformality::Vec2> fd_e1;
  std::vector<conformality::Vec2> fd_e2;
};

#endif  // CONFORMALITY_TRACERS_HPP_
