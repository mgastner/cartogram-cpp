#ifndef REPAIR_RINGS_HPP_
#define REPAIR_RINGS_HPP_

#include "geo_div.hpp"
#include <vector>

// Whether a ring is simple, tested with exact predicates. The repair below and
// InsetState::is_simple() decide with this test.
bool is_simple_with_exact_predicates(const Polygon &ring);

// Repairs the rings of geo_divs where rounding has made a ring cross or touch
// itself (see src/inset_state/repair_rings.cpp). The rings are numbered in
// the order of the GeoDivs and their polygons with holes, each outer boundary
// before its holes. Both functions return whether every ring is simple.

// After densification: geo_divs holds the densified rings and before the
// rings before densification, and inserted[r][k] tells whether
// densification inserted point k of ring r. Only rings that were simple
// before are repaired, and only inserted points are moved or removed.
bool repair_rings_after_densification(
  std::vector<GeoDiv> &geo_divs,
  const std::vector<GeoDiv> &before,
  std::vector<std::vector<bool>> &inserted);

// After the projection: any point may be moved, and none is removed
bool repair_rings_after_projection(std::vector<GeoDiv> &geo_divs);

#endif  // REPAIR_RINGS_HPP_
