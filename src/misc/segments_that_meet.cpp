#include "segments_that_meet.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <algorithm>
#include <numeric>

// Kernel with exact predicates. Its points have the same double coordinates
// as the project's points, so only the predicates differ.
using Epick = CGAL::Exact_predicates_inexact_constructions_kernel;

static std::vector<Epick::Point_2> with_exact_predicates(const Polygon &ring)
{
  std::vector<Epick::Point_2> pts;
  pts.reserve(ring.size());
  for (const Point &p : ring)
    pts.emplace_back(p.x(), p.y());
  return pts;
}

bool is_simple_with_exact_predicates(const Polygon &ring)
{
  const std::vector<Epick::Point_2> pts = with_exact_predicates(ring);
  return CGAL::is_simple_2(pts.begin(), pts.end(), Epick());
}

// Whether the segments (a, s) and (s, b), which are neighbors in a ring,
// overlap beyond their common end s
static bool fold_back(
  const Epick::Point_2 &a,
  const Epick::Point_2 &s,
  const Epick::Point_2 &b)
{
  return CGAL::collinear(a, s, b) &&
         !CGAL::collinear_are_strictly_ordered_along_line(a, s, b);
}

// The segments are swept in the order of their smallest x-coordinates
std::vector<SegmentsThatMeet> segments_that_meet(const Polygon &ring)
{
  std::vector<SegmentsThatMeet> found;
  const std::size_t n = ring.size();
  if (n < 3)
    return found;
  const std::vector<Epick::Point_2> pts = with_exact_predicates(ring);
  auto next = [n](std::size_t k) {
    return (k + 1 == n) ? 0 : k + 1;
  };
  auto lo_x = [&](std::size_t k) {
    return std::min(pts[k].x(), pts[next(k)].x());
  };
  auto hi_x = [&](std::size_t k) {
    return std::max(pts[k].x(), pts[next(k)].x());
  };
  auto lo_y = [&](std::size_t k) {
    return std::min(pts[k].y(), pts[next(k)].y());
  };
  auto hi_y = [&](std::size_t k) {
    return std::max(pts[k].y(), pts[next(k)].y());
  };
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    return lo_x(a) < lo_x(b);
  });
  for (std::size_t a = 0; a < n; ++a) {
    const std::size_t i = order[a];
    for (std::size_t b = a + 1; b < n && !(lo_x(order[b]) > hi_x(i)); ++b) {
      const std::size_t j = order[b];
      if (lo_y(j) > hi_y(i) || lo_y(i) > hi_y(j))
        continue;
      SegmentsThatMeet meet;
      meet.i = i;
      meet.j = j;
      if (next(i) == j || next(j) == i) {
        meet.fold_back = (next(i) == j)
                           ? fold_back(pts[i], pts[j], pts[next(j)])
                           : fold_back(pts[j], pts[i], pts[next(i)]);
        if (meet.fold_back)
          found.push_back(meet);
        continue;
      }
      if (!CGAL::do_intersect(
            Epick::Segment_2(pts[i], pts[next(i)]),
            Epick::Segment_2(pts[j], pts[next(j)])))
        continue;

      // Segments that are not neighbors in the ring and have exactly one end
      // in common, without overlapping beyond it, touch only there
      std::size_t n_common = 0;
      for (const std::size_t ku : {i, next(i)}) {
        for (const std::size_t kv : {j, next(j)}) {
          if (pts[ku] == pts[kv]) {
            ++n_common;
            meet.u = ku;
            meet.v = kv;
          }
        }
      }
      meet.touch_at_common_end =
        n_common == 1 && !fold_back(
                           pts[meet.u == i ? next(i) : i],
                           pts[meet.u],
                           pts[meet.v == j ? next(j) : j]);
      found.push_back(meet);
    }
  }
  return found;
}
