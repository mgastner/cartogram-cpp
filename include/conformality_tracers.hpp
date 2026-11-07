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
  // x + fd_delta e_i (fd_e1, fd_e2) and at x + (fd_delta / 2) e_i (fd_e1h,
  // fd_e2h). Forward differences are exact for the piecewise-affine T0
  // unless a pair straddles a piece; after the smooth correction A, the
  // Richardson combination of both step sizes cancels the O(delta) error.
  double fd_delta{1e-4};
  std::vector<std::size_t> fd_parent;
  std::vector<conformality::Vec2> fd_e1;
  std::vector<conformality::Vec2> fd_e2;
  std::vector<conformality::Vec2> fd_e1h;
  std::vector<conformality::Vec2> fd_e2h;

  // Square graticule for figures: lines of constant x and y in the initial
  // lattice, every grid_spacing lattice units, sampled at
  // grid_samples_per_unit points per lattice unit. grid_lines holds the
  // current positions, grid_lines_t0 a snapshot after the last pass.
  double grid_spacing{8.0};
  unsigned int grid_samples_per_unit{4};
  std::vector<std::vector<conformality::Vec2>> grid_lines;
  std::vector<std::vector<conformality::Vec2>> grid_lines_t0;

  // Distortion field for figures, on cells of field_cell x field_cell units
  // of the initial lattice (field_nx x field_ny cells, the last row and
  // column possibly narrower). field_corners holds the current images of
  // the cell corners, row by row from y = 0, and field_corners_t0 a
  // snapshot after the last pass. field_e_t0 and field_m_t0 hold the mean
  // ln(a/b) and mean excess MIPS energy of T0 over the evaluation nodes in
  // each cell, in the same order.
  double field_cell{2.0};
  unsigned int field_nx{0};
  unsigned int field_ny{0};
  std::vector<conformality::Vec2> field_corners;
  std::vector<conformality::Vec2> field_corners_t0;
  std::vector<double> field_e_t0;
  std::vector<double> field_m_t0;

  // Rings of the input polygons at seeding time, with the index of their
  // GeoDiv, for figures of the original map
  std::vector<std::vector<conformality::Vec2>> input_rings;
  std::vector<int> input_ring_region;

  // After the optimizer ran: the output polygons of T0, densified, with the
  // correction A applied (for figures of T = A o T0)
  bool optimized{false};
  std::vector<std::vector<conformality::Vec2>> t_rings;
  std::vector<int> t_ring_region;
};

#endif  // CONFORMALITY_TRACERS_HPP_
