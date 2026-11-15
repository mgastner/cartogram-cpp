#ifndef SEGMENTS_THAT_MEET_HPP_
#define SEGMENTS_THAT_MEET_HPP_

#include "cgal_typedef.hpp"
#include <cstddef>
#include <vector>

// Two segments of a ring that meet, other than neighbors at their common end.
// Segment k runs from ring[k] to the next point of the ring.
struct SegmentsThatMeet {
  std::size_t i = 0;
  std::size_t j = 0;

  // Whether the segments are neighbors in the ring and overlap beyond their
  // common end
  bool fold_back = false;

  // Whether the segments are not neighbors, have exactly one end in common,
  // and do not overlap beyond it, so that they touch only there. The common
  // end occurs twice in the ring, as ring[u] and ring[v].
  bool touch_at_common_end = false;
  std::size_t u = 0;
  std::size_t v = 0;
};

// The pairs of segments of the ring that meet, found with exact predicates
std::vector<SegmentsThatMeet> segments_that_meet(const Polygon &ring);

// Whether the ring is simple according to exact predicates. The simplicity
// test of the Polygon type rounds.
bool is_simple_with_exact_predicates(const Polygon &ring);

#endif  // SEGMENTS_THAT_MEET_HPP_
