#ifndef REPAIR_RINGS_HPP_
#define REPAIR_RINGS_HPP_

#include "geo_div.hpp"
#include <vector>

// Repairs the rings of geo_divs where rounding has made a ring cross or touch
// itself (see src/inset_state/repair_rings.cpp). The rings are numbered in
// the order of the GeoDivs and their polygons with holes, each outer boundary
// before its holes. Both functions return whether every ring passes the
// simplicity test of the Polygon type, which is also the test of
// InsetState::is_simple().

// After densification: geo_divs holds the densified rings and before the
// rings before densification, and inserted[r][k] tells whether
// densification inserted point k of ring r. Only rings that passed the test
// before are repaired, and only inserted points are moved or removed.
bool repair_rings_after_densification(
  std::vector<GeoDiv> &geo_divs,
  const std::vector<GeoDiv> &before,
  std::vector<std::vector<bool>> &inserted);

// After the projection: any point may be moved, and none is removed
bool repair_rings_after_projection(std::vector<GeoDiv> &geo_divs);

#endif  // REPAIR_RINGS_HPP_
