#include "repair_rings.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

// Kernel with exact predicates for checking simplicity. Its points have the
// same double coordinates as the project's points, so only the predicates
// differ.
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

// Whether x lies strictly between s and t on the segment from s to t. The
// test uses exact predicates.
static bool strictly_inside(const Point &s, const Point &x, const Point &t)
{
  const Epick::Point_2 es(s.x(), s.y());
  const Epick::Point_2 ex(x.x(), x.y());
  const Epick::Point_2 et(t.x(), t.y());
  return CGAL::collinear(es, ex, et) &&
         CGAL::collinear_are_strictly_ordered_along_line(es, ex, et);
}

// Whether the segments from s to a and from s to b, which are segments of two
// rings at a point s that both rings contain, overlap beyond s, as on a
// shared boundary. The test uses exact predicates.
static bool segments_overlap(const Point &a, const Point &s, const Point &b)
{
  return a != s && b != s &&
         fold_back(
           Epick::Point_2(a.x(), a.y()),
           Epick::Point_2(s.x(), s.y()),
           Epick::Point_2(b.x(), b.y()));
}

// The pairs of segments of a ring that meet, other than neighbors at their
// common end, and the movable points at their ends. Segment k of a ring runs
// from point k to the next point.
struct SegmentsThatMeet {
  std::size_t pairs = 0;

  // Sorted indices k with movable[k]
  std::vector<std::size_t> ends;
};

// Adds the segments i and j of the ring to found if they meet, which is
// tested with exact predicates. If the two segments only touch at an end that
// occurs twice in the ring, only that end counts as an end of the pair.
static void add_if_segments_meet(
  const Polygon &ring,
  const std::vector<bool> &movable,
  const std::size_t i,
  const std::size_t j,
  SegmentsThatMeet &found)
{
  const std::size_t n = ring.size();
  auto next = [n](std::size_t k) {
    return (k + 1 == n) ? 0 : k + 1;
  };
  auto pt = [&ring](std::size_t k) {
    return Epick::Point_2(ring[k].x(), ring[k].y());
  };
  bool meet;
  if (next(i) == j) {
    meet = fold_back(pt(i), pt(j), pt(next(j)));
  } else if (next(j) == i) {
    meet = fold_back(pt(j), pt(i), pt(next(i)));
  } else {
    meet = CGAL::do_intersect(
      Epick::Segment_2(pt(i), pt(next(i))),
      Epick::Segment_2(pt(j), pt(next(j))));
  }
  if (!meet)
    return;
  ++found.pairs;

  // Segments that are not neighbors in the ring and have exactly one end in
  // common, without overlapping beyond it, touch only there
  std::size_t n_common = 0, u = 0, v = 0;
  if (next(i) != j && next(j) != i) {
    for (const std::size_t ku : {i, next(i)}) {
      for (const std::size_t kv : {j, next(j)}) {
        if (ring[ku] == ring[kv]) {
          ++n_common;
          u = ku;
          v = kv;
        }
      }
    }
  }
  const bool touch_at_common_end =
    n_common == 1 &&
    !fold_back(pt(u == i ? next(i) : i), pt(u), pt(v == j ? next(j) : j));
  if (touch_at_common_end) {
    for (const std::size_t k : {u, v}) {
      if (movable[k])
        found.ends.push_back(k);
    }
    return;
  }
  for (const std::size_t k : {i, next(i), j, next(j)}) {
    if (movable[k])
      found.ends.push_back(k);
  }
}

// Bounding box of segment k of a ring
struct SegmentBox {
  double lo_x, hi_x, lo_y, hi_y;
};

static SegmentBox box_of(const Polygon &ring, const std::size_t k)
{
  const Point &a = ring[k];
  const Point &b = ring[(k + 1 == ring.size()) ? 0 : k + 1];
  return {
    std::min(a.x(), b.x()),
    std::max(a.x(), b.x()),
    std::min(a.y(), b.y()),
    std::max(a.y(), b.y())};
}

static void sort_ends(SegmentsThatMeet &found)
{
  std::sort(found.ends.begin(), found.ends.end());
  found.ends.erase(
    std::unique(found.ends.begin(), found.ends.end()),
    found.ends.end());
}

// The pairs of segments of the ring that meet. The segments are swept in the
// order of their smallest x-coordinates.
static SegmentsThatMeet segments_that_meet(
  const Polygon &ring,
  const std::vector<bool> &movable)
{
  SegmentsThatMeet found;
  const std::size_t n = ring.size();
  if (n < 3)
    return found;
  std::vector<SegmentBox> boxes;
  boxes.reserve(n);
  for (std::size_t k = 0; k < n; ++k)
    boxes.push_back(box_of(ring, k));
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    return boxes[a].lo_x < boxes[b].lo_x;
  });
  for (std::size_t a = 0; a < n; ++a) {
    const std::size_t i = order[a];
    for (std::size_t b = a + 1;
         b < n && !(boxes[order[b]].lo_x > boxes[i].hi_x);
         ++b) {
      const std::size_t j = order[b];
      if (boxes[j].lo_y > boxes[i].hi_y || boxes[i].lo_y > boxes[j].hi_y)
        continue;
      add_if_segments_meet(ring, movable, i, j, found);
    }
  }
  sort_ends(found);
  return found;
}

// The pairs of segments of the ring that meet among the pairs with a segment
// at point k, which are the only pairs that moving point k can change
static SegmentsThatMeet segments_that_meet_at(
  const Polygon &ring,
  const std::vector<bool> &movable,
  const std::size_t k)
{
  SegmentsThatMeet found;
  const std::size_t n = ring.size();
  if (n < 3)
    return found;
  const std::size_t prev = (k + n - 1) % n;
  for (const std::size_t i : {prev, k}) {
    const SegmentBox bi = box_of(ring, i);
    for (std::size_t j = 0; j < n; ++j) {
      if (j == i || (i == k && j == prev))
        continue;
      const SegmentBox bj = box_of(ring, j);
      if (
        bj.lo_x > bi.hi_x || bi.lo_x > bj.hi_x || bj.lo_y > bi.hi_y ||
        bi.lo_y > bj.hi_y)
        continue;
      add_if_segments_meet(ring, movable, i, j, found);
    }
  }
  sort_ends(found);
  return found;
}

// The eight doubles next to p, one unit in the last place away in x, in y,
// or in both. Those that differ from p in one coordinate come first.
static std::array<Point, 8> neighboring_points(const Point &p)
{
  const double inf = std::numeric_limits<double>::infinity();
  const double xl = std::nextafter(p.x(), -inf);
  const double xh = std::nextafter(p.x(), inf);
  const double yl = std::nextafter(p.y(), -inf);
  const double yh = std::nextafter(p.y(), inf);
  return {
    Point(p.x(), yl),
    Point(p.x(), yh),
    Point(xl, p.y()),
    Point(xh, p.y()),
    Point(xl, yl),
    Point(xh, yl),
    Point(xl, yh),
    Point(xh, yh)};
}

// Rounding can make a ring cross or touch itself where another part of the
// ring is closer than a unit in the last place of the coordinates, for
// example the other side of a needle one unit in the last place wide:
// - Densification inserts the crossing of a polygon edge with a triangle
//   edge, rounded to doubles, which can lie up to half a unit in the last
//   place off the polygon edge: across the other side of the needle, or on a
//   vertex of that side.
// - The projection maps each point with the affine map of its triangle and
//   rounds the image, which can put the image of a point across or onto the
//   image of a nearby segment, although no triangle is inverted.
// Every ring is checked for simplicity with exact predicates, as in
// InsetState::is_simple(). In a ring that fails, the movable points at the
// ends of the segments that meet are found. Each such point is moved to the
// one of the eight doubles next to it that leaves the fewest pairs of
// segments that meet, if that is fewer than before the move.
//
// A point moves together with its copies: the occurrences of the same point
// in other rings with a segment that overlaps a segment at the point, as on a
// shared boundary, also where one ring has vertices on that boundary that
// the other lacks. A move must keep every such pair of segments overlapping,
// so that shared boundaries stay shared; a neighbor that would separate them
// is not a candidate. Either every copy moves, or none does.
//
// Moving a point by a unit in the last place keeps each densified edge within
// about a unit in the last place of a triangle, so the projection still maps
// the polygon edge piecewise affinely and follows the bends of its exact
// image. Removing an inserted point would replace the image of the polygon
// edge between the neighbors of the point by a straight segment, which can
// cross the other side of a thin part of a polygon. After densification, an
// inserted point is therefore removed only if no move reduces the number of
// pairs of segments that meet, neither of the point nor of any other point in
// the same round. Each ring keeps its vertices from before densification, so
// the worst case for a ring is its state before densification. After the
// projection, no point is removed, and is_simple() reports a ring that the
// moves do not repair.
//
// The rings that changed are checked again until no ring changes. Each move
// reduces the number of pairs of segments that meet, and each point is
// removed at most once, so the repair ends.
static bool repair_rings(
  std::vector<GeoDiv> &geo_divs,
  const std::vector<GeoDiv> *before,
  std::vector<std::vector<bool>> &movable,
  const std::string &stage)
{
  std::vector<Polygon *> rings;
  for (auto &gd : geo_divs) {
    for (auto &pwh : gd.ref_to_polygons_with_holes()) {
      rings.push_back(&pwh.outer_boundary());
      for (auto &h : pwh.holes())
        rings.push_back(&h);
    }
  }
  std::vector<const Polygon *> rings_before;
  if (before) {
    for (const auto &gd : *before) {
      for (const auto &pwh : gd.polygons_with_holes()) {
        rings_before.push_back(&pwh.outer_boundary());
        for (const auto &h : pwh.holes())
          rings_before.push_back(&h);
      }
    }
    assert(rings_before.size() == rings.size());
  }
  assert(movable.size() == rings.size());

  // A point of a ring and the index of the ring
  using Occurrence = std::pair<std::size_t, std::size_t>;
  auto point_at = [&](const Occurrence &o) -> Point & {
    return (*rings[o.first])[o.second];
  };
  auto neighbors_of = [&](const Occurrence &o) {
    const Polygon &ring = *rings[o.first];
    const std::size_t n = ring.size();
    return std::array<Point, 2>{
      ring[(o.second + n - 1) % n],
      ring[(o.second + 1) % n]};
  };

  std::vector<bool> fails(rings.size(), false);
  std::vector<std::size_t> to_check(rings.size());
  std::iota(to_check.begin(), to_check.end(), std::size_t{0});
  std::size_t n_moved = 0;
  std::unordered_set<Point> removed;
  while (!to_check.empty()) {
    std::vector<Occurrence> found;
    for (const std::size_t r : to_check) {
      const Polygon &ring = *rings[r];
      fails[r] = !is_simple_with_exact_predicates(ring);
      if (
        !fails[r] ||
        (before && !is_simple_with_exact_predicates(*rings_before[r])))
        continue;
      for (const std::size_t k : segments_that_meet(ring, movable[r]).ends)
        found.emplace_back(r, k);
    }
    if (found.empty())
      break;

    // Where the points found occur in movable positions of any ring
    std::unordered_set<Point> points;
    for (const Occurrence &o : found)
      points.insert(point_at(o));
    std::unordered_map<Point, std::vector<Occurrence>> at;
    for (std::size_t r = 0; r < rings.size(); ++r) {
      for (std::size_t k = 0; k < rings[r]->size(); ++k) {
        if (movable[r][k] && points.count((*rings[r])[k]) > 0)
          at[(*rings[r])[k]].emplace_back(r, k);
      }
    }

    // An occurrence and its copies: the occurrences of the same point in
    // other rings with a segment that overlaps a segment at the occurrence or
    // at another copy. An occurrence that a move in this round has moved is no
    // longer a copy.
    auto with_copies = [&](const Occurrence &o) {
      const Point p = point_at(o);
      std::vector<Occurrence> group = {o};
      for (std::size_t g = 0; g < group.size(); ++g) {
        const Occurrence m = group[g];
        for (const Occurrence &q : at[p]) {
          if (
            q.first == m.first || point_at(q) != p ||
            std::find(group.begin(), group.end(), q) != group.end())
            continue;
          bool linked = false;
          for (const Point &a : neighbors_of(m)) {
            for (const Point &b : neighbors_of(q))
              linked = linked || segments_overlap(a, p, b);
          }
          if (linked)
            group.push_back(q);
        }
      }
      return group;
    };

    // Where the boundaries of other rings run along the segments at the
    // occurrences of a group: the vertices of other rings that lie strictly
    // inside a segment at an occurrence, and the segments of other rings that
    // have the point strictly inside. A ring that contains the point as a
    // vertex and also has it strictly inside a segment touches itself there,
    // which is what the repair removes, so such segments do not count. Each
    // contact is a triple (s, x, t) with x strictly inside the segment from s
    // to t, in which s, x, or t is the point. A move to c keeps shared
    // boundaries shared, also where one ring has vertices on them that the
    // other lacks, only if every triple still holds with c in place of the
    // point.
    auto contacts_of = [&](const std::vector<Occurrence> &group) {
      const Point p = point_at(group.front());
      std::vector<std::array<Point, 3>> contacts;
      for (const Occurrence &m : group) {
        const auto nb = neighbors_of(m);
        for (std::size_t r = 0; r < rings.size(); ++r) {
          if (r == m.first)
            continue;
          const Polygon &ring = *rings[r];
          bool contains_p = false;
          std::vector<std::array<Point, 3>> segments_through_p;
          for (std::size_t k = 0; k < ring.size(); ++k) {
            const Point &v = ring[k];
            const Point &w = ring[(k + 1 == ring.size()) ? 0 : k + 1];
            contains_p = contains_p || v == p;
            for (const Point &a : nb) {
              if (v != p && strictly_inside(a, v, p))
                contacts.push_back({a, v, p});
            }
            if (strictly_inside(v, p, w))
              segments_through_p.push_back({v, p, w});
          }
          if (!contains_p) {
            contacts.insert(
              contacts.end(),
              segments_through_p.begin(),
              segments_through_p.end());
          }
        }
      }
      return contacts;
    };

    // Number of pairs of segments that meet in the rings of a group, up to
    // pairs that a move of the group does not change: in a ring with one
    // occurrence, only the pairs with a segment at the occurrence count
    auto pairs_in = [&](const std::vector<Occurrence> &group) {
      std::map<std::size_t, std::vector<std::size_t>> by_ring;
      for (const Occurrence &q : group)
        by_ring[q.first].push_back(q.second);
      std::size_t total = 0;
      for (const auto &[r, ks] : by_ring) {
        total += (ks.size() == 1)
                   ? segments_that_meet_at(*rings[r], movable[r], ks[0]).pairs
                   : segments_that_meet(*rings[r], movable[r]).pairs;
      }
      return total;
    };

    // The groups that no move repairs. After densification, they are removed
    // if no group moves in this round; otherwise, the rings are checked again.
    std::vector<Occurrence> done;
    std::vector<Occurrence> stuck;
    std::size_t n_moved_now = 0;
    std::vector<bool> changed_now(rings.size(), false);
    for (const Occurrence &o : found) {
      if (std::find(done.begin(), done.end(), o) != done.end())
        continue;

      // A move earlier in this round may have repaired the occurrence
      const std::vector<std::size_t> ks =
        segments_that_meet_at(*rings[o.first], movable[o.first], o.second)
          .ends;
      if (!std::binary_search(ks.begin(), ks.end(), o.second))
        continue;
      const std::vector<Occurrence> group = with_copies(o);
      done.insert(done.end(), group.begin(), group.end());
      const Point p = point_at(o);
      const std::vector<std::array<Point, 3>> contacts = contacts_of(group);

      // Move the group to the neighbor that leaves the fewest pairs of
      // segments that meet, if that is fewer than now, among the neighbors
      // that keep shared boundaries shared. Either every occurrence in the
      // group moves, or none does.
      std::size_t fewest = pairs_in(group);
      Point best = p;
      for (const Point &c : neighboring_points(p)) {
        const bool keeps_contacts = std::all_of(
          contacts.begin(),
          contacts.end(),
          [&p, &c](std::array<Point, 3> stx) {
            for (Point &e : stx) {
              if (e == p)
                e = c;
            }
            return strictly_inside(stx[0], stx[1], stx[2]);
          });
        if (!keeps_contacts)
          continue;
        for (const Occurrence &q : group)
          point_at(q) = c;
        const std::size_t n_pairs = pairs_in(group);
        if (n_pairs < fewest) {
          fewest = n_pairs;
          best = c;
        }
      }
      for (const Occurrence &q : group)
        point_at(q) = best;
      if (best != p) {
        ++n_moved_now;
        for (const Occurrence &q : group)
          changed_now[q.first] = true;
      } else {
        stuck.insert(stuck.end(), group.begin(), group.end());
      }
    }
    n_moved += n_moved_now;
    std::vector<Occurrence> to_remove;
    if (n_moved_now == 0 && before)
      to_remove = stuck;

    // Remove the points that could not be moved, with their copies
    std::vector<std::vector<std::size_t>> remove_at(rings.size());
    for (const Occurrence &o : to_remove)
      remove_at[o.first].push_back(o.second);
    for (std::size_t r = 0; r < rings.size(); ++r) {
      if (remove_at[r].empty())
        continue;
      std::sort(remove_at[r].begin(), remove_at[r].end());
      Polygon &ring = *rings[r];
      Polygon kept;
      std::vector<bool> kept_movable;
      for (std::size_t k = 0; k < ring.size(); ++k) {
        if (std::binary_search(remove_at[r].begin(), remove_at[r].end(), k)) {
          removed.insert(ring[k]);
        } else {
          kept.push_back(ring[k]);
          kept_movable.push_back(movable[r][k]);
        }
      }
      ring = std::move(kept);
      movable[r] = std::move(kept_movable);
      changed_now[r] = true;
    }

    to_check.clear();
    for (std::size_t r = 0; r < rings.size(); ++r) {
      if (changed_now[r])
        to_check.push_back(r);
    }
  }
  if (n_moved > 0) {
    std::cerr << "Moved " << n_moved
              << " point(s) by a unit in the last place where rounding in "
              << stage << " made a ring cross itself" << std::endl;
  }
  if (!removed.empty()) {
    std::cerr << "Removed " << removed.size()
              << " densification point(s) that made a ring cross itself"
              << std::endl;
  }
  return std::none_of(fails.begin(), fails.end(), [](bool f) {
    return f;
  });
}

bool repair_rings_after_densification(
  std::vector<GeoDiv> &geo_divs,
  const std::vector<GeoDiv> &before,
  std::vector<std::vector<bool>> &inserted)
{
  return repair_rings(geo_divs, &before, inserted, "densification");
}

bool repair_rings_after_projection(std::vector<GeoDiv> &geo_divs)
{
  std::vector<std::vector<bool>> movable;
  for (const auto &gd : geo_divs) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      movable.emplace_back(pwh.outer_boundary().size(), true);
      for (const auto &h : pwh.holes())
        movable.emplace_back(h.size(), true);
    }
  }
  return repair_rings(geo_divs, nullptr, movable, "the projection");
}
