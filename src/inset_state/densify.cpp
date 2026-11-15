#include "inset_state.hpp"
#include "round_point.hpp"
#include <CGAL/AABB_segment_primitive_2.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

using AABB_Prim =
  CGAL::AABB_segment_primitive_2<Scd, std::vector<Segment>::const_iterator>;
using AABB_Traits = CGAL::AABB_traits_2<Scd, AABB_Prim>;
using AABB_Tree = CGAL::AABB_tree<AABB_Traits>;

using SegCache = std::unordered_map<Segment, std::vector<Point>>;

// Edge classification of AABB tree segments. AABB tree primitives give the
// index of the segment in the vector that we passed to create the tree, so
// we can classify edges by their index in that vector. A crossing with a
// horizontal or vertical edge is placed exactly on its lattice line.
//
// H    (0): horizontal edge, y1 == y2
// V    (1): vertical edge,   x1 == x2
// OTHER(2): all other edges (diagonals of slope +1 or -1, and edges of
//           other slopes where quadtree leaves of different sizes meet)
enum class EdgeKind : uint8_t { H = 0, V = 1, OTHER = 2 };

static inline EdgeKind classify_edge(
  const uint32_t x0,
  const uint32_t y0,
  const uint32_t x1,
  const uint32_t y1)
{
  if (y0 == y1)
    return EdgeKind::H;
  if (x0 == x1)
    return EdgeKind::V;
  return EdgeKind::OTHER;
}

// Error-free transformations: s + e == a + b (Knuth's TwoSum) and
// p + e == a * b (TwoProduct with a fused multiply-add), both exactly. They
// need round-to-nearest arithmetic without contraction or reassociation; the
// project compiles with -ffp-contract=off and without -ffast-math.
static void two_sum(const double a, const double b, double &s, double &e)
{
  s = a + b;
  const double b_virtual = s - a;
  e = (a - (s - b_virtual)) + (b - b_virtual);
}

static void two_product(const double a, const double b, double &p, double &e)
{
  p = a * b;
  e = std::fma(a, b, -p);
}

// Signed offset nx * px + ny * py - c of the point (px, py) from the line
// nx * x + ny * y = c. The sum is evaluated as if in twice the working
// precision and then rounded (algorithm Dot2 of Ogita, Rump, and Oishi,
// 2005), so the offset keeps its relative accuracy even when the point lies
// within a rounding error of the line. Evaluated directly, the sum cancels:
// for a point 2e-12 from the line x + y = 32, rounding px + py alone
// changes the offset by about 0.1%.
static double offset_from_line(
  const double nx,
  const double ny,
  const double c,
  const double px,
  const double py)
{
  double p1, e1, p2, e2, s1, e3, s2, e4;
  two_product(nx, px, p1, e1);
  two_product(ny, py, p2, e2);
  two_sum(p1, p2, s1, e3);
  two_sum(s1, -c, s2, e4);
  return s2 + (e1 + e2 + e3 + e4);
}

static std::vector<Point> densification_points_with_edge_tree(
  const Point &pt1,
  const Point &pt2,
  const AABB_Tree &tree,
  const std::vector<Segment> &edge_vec,
  const std::vector<EdgeKind> &edge_kind,
  SegCache *cache)
{
  // Check in cache first. We have lots of shared boundaries, so caching will
  // be very useful. Basically we will be doing almost half of the work we
  // would normally do
  if (cache) {
    auto it = cache->find({pt1, pt2});
    if (it != cache->end())
      return it->second;

    // Hole and outer boundaries have reverse vertex order in the polygon, so
    // if we processed (pt1, pt2), the shared boundary of the hole of some
    // other GeoDiv will be (pt2, pt1). So, we check both
    auto ir = cache->find({pt2, pt1});
    if (ir != cache->end()) {
      std::vector<Point> rev = ir->second;
      std::reverse(rev.begin(), rev.end());
      return rev;
    }
  }

  const Segment query(pt1, pt2);
  const double ax = pt1.x(), ay = pt1.y();
  const double bx = pt2.x(), by = pt2.y();
  const double dax = bx - ax, day = by - ay;

  // Collect primitives that intersect the query (not just bbox hits)
  std::vector<AABB_Tree::Primitive_id> prims;
  prims.reserve(8);
  tree.all_intersected_primitives(query, std::back_inserter(prims));

  // Always include endpoints. The parameter t along [pt1, pt2] orders the
  // points.
  std::vector<std::pair<double, Point>> hits;
  hits.reserve(prims.size() + 2);
  hits.emplace_back(0.0, pt1);
  hits.emplace_back(1.0, pt2);

  auto within_closed = [](double v, double a, double b) {
    const double lo = std::min(a, b);
    const double hi = std::max(a, b);
    return less_than_equal(lo, v) && less_than_equal(v, hi);
  };

  for (const auto &pid : prims) {

    // Primitive_id is an iterator into edge_vec
    const std::ptrdiff_t idx = pid - edge_vec.begin();

    assert(idx >= 0 && static_cast<size_t>(idx) < edge_vec.size());

    const Segment &E = *pid;
    const EdgeKind kind = edge_kind[static_cast<size_t>(idx)];

    const double ex0 = E.source().x(), ey0 = E.source().y();
    const double ex1 = E.target().x(), ey1 = E.target().y();

    // The triangle edge lies on the line nx * x + ny * y = c. Its ends are
    // corners of quadtree leaves, which have integer coordinates, so nx, ny,
    // and c are integers and exact.
    const double nx = ey0 - ey1;
    const double ny = ex1 - ex0;
    const double c = nx * ex0 + ny * ey0;
    const double fa = offset_from_line(nx, ny, c, ax, ay);
    const double fb = offset_from_line(nx, ny, c, bx, by);

    // If the polygon edge lies within dbl_resolution of the line along its
    // whole length, the position of the crossing along the edge is
    // ill-conditioned, and rounding the coordinates of an inserted point
    // could move the point across the other side of a thin part of the
    // polygon. Such an edge gets no point from this triangle edge. The
    // projected edge then differs from the exact image of the edge by at
    // most dbl_resolution times the jump of the Jacobian across the line.
    const double tolerance = dbl_resolution * std::sqrt(nx * nx + ny * ny);
    if (std::fabs(fa) <= tolerance && std::fabs(fb) <= tolerance)
      continue;

    // The ends must lie strictly on opposite sides. An end on the line is
    // already among the hits.
    if (!((fa < 0.0 && fb > 0.0) || (fa > 0.0 && fb < 0.0)))
      continue;

    // Rounding keeps the ratio of two offsets of opposite signs within
    // [0, 1]. A ratio rounded to 0 or 1 is an end of the edge.
    const double t = fa / (fa - fb);
    if (!(t > 0.0 && t < 1.0))
      continue;

    // A crossing with a horizontal or vertical edge lies exactly on its
    // lattice line
    const double x = (kind == EdgeKind::V) ? ex0 : ax + t * dax;
    const double y = (kind == EdgeKind::H) ? ey0 : ay + t * day;
    if (within_closed(x, ex0, ex1) && within_closed(y, ey0, ey1))
      hits.emplace_back(t, Point(x, y));
  }

  // Sort by t; clamp and unique
  std::sort(hits.begin(), hits.end(), [](auto const &a, auto const &b) {
    return a.first < b.first;
  });

  std::vector<Point> out_unique;
  out_unique.reserve(hits.size());
  for (const auto &hit : hits) {
    out_unique.push_back(hit.second);
  }

  out_unique.erase(
    std::unique(out_unique.begin(), out_unique.end()),
    out_unique.end());

  if (cache) {
    (*cache)[{pt1, pt2}] = out_unique;
  }
  return out_unique;
}

// Appends the points of one densified polygon edge, without its last point,
// to a ring, and records for each point whether densification inserted it
static void append_densified_edge(
  Polygon &ring,
  std::vector<bool> &inserted,
  const std::vector<Point> &seg_pts)
{
  if (seg_pts.size() < 2)
    return;
  ring.insert(ring.end(), seg_pts.begin(), seg_pts.end() - 1);
  inserted.push_back(false);
  inserted.insert(inserted.end(), seg_pts.size() - 2, true);
}

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

static bool is_simple_with_exact_predicates(const Polygon &ring)
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

// Collects the points ring[k] with inserted[k] that are ends of a segment of
// the ring that meets another segment of the ring, other than its neighbors
// at their common end. If the two segments only touch at an end that occurs
// twice in the ring, only that end is collected. The test uses exact
// predicates. The segments are swept in the order of their smallest
// x-coordinates.
static void densification_points_on_intersecting_segments(
  const Polygon &ring,
  const std::vector<bool> &inserted,
  std::unordered_set<Point> &found)
{
  const std::size_t n = ring.size();
  if (n < 3)
    return;
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
      bool meet;
      if (next(i) == j) {
        meet = fold_back(pts[i], pts[j], pts[next(j)]);
      } else if (next(j) == i) {
        meet = fold_back(pts[j], pts[i], pts[next(i)]);
      } else {
        meet = CGAL::do_intersect(
          Epick::Segment_2(pts[i], pts[next(i)]),
          Epick::Segment_2(pts[j], pts[next(j)]));
      }
      if (!meet)
        continue;

      // Segments that are not neighbors in the ring and have exactly one
      // end in common, without overlapping beyond it, touch only there
      std::size_t n_common = 0, u = 0, v = 0;
      if (next(i) != j && next(j) != i) {
        for (const std::size_t ku : {i, next(i)}) {
          for (const std::size_t kv : {j, next(j)}) {
            if (pts[ku] == pts[kv]) {
              ++n_common;
              u = ku;
              v = kv;
            }
          }
        }
      }
      const bool touch_at_common_end =
        n_common == 1 && !fold_back(
                           pts[u == i ? next(i) : i],
                           pts[u],
                           pts[v == j ? next(j) : j]);
      if (touch_at_common_end) {
        for (const std::size_t k : {u, v}) {
          if (inserted[k])
            found.insert(ring[k]);
        }
        continue;
      }
      for (const std::size_t k : {i, next(i), j, next(j)}) {
        if (inserted[k])
          found.insert(ring[k]);
      }
    }
  }
}

// Rounding the coordinates of an inserted point moves the point off its
// polygon edge by up to half a unit in the last place. Where another part of
// the ring is that close, for example the other side of a needle one unit in
// the last place wide, the ring can cross itself, or both sides can receive
// the same point, so that the ring touches itself. Every ring is therefore
// checked with the simplicity test of the Polygon type, which is also the
// test of InsetState::is_simple(). If a ring that passed the test before
// densification fails it now, the inserted points at the ends of the
// segments that meet, found with exact predicates, are removed from every
// ring wherever densification inserted them (inserted[r][k] tells whether
// it inserted point k of ring r), so that shared boundaries stay identical.
// If no such segments are found, all inserted points of the ring are
// removed. The rings that lost points are checked again, with the test of
// the Polygon type and with exact predicates, until no ring changes. Each
// ring keeps its vertices from before densification, so the worst case for
// a ring is its state before densification. Returns whether every ring
// passes the test of the Polygon type.
static bool drop_densification_points_that_break_simplicity(
  std::vector<GeoDiv> &geo_divs_dens,
  const std::vector<GeoDiv> &geo_divs,
  std::vector<std::vector<bool>> &inserted)
{
  // Each densified ring with the ring before densification
  std::vector<std::pair<Polygon *, const Polygon *>> rings;
  for (std::size_t g = 0; g < geo_divs.size(); ++g) {
    auto &pwhs_dens = geo_divs_dens[g].ref_to_polygons_with_holes();
    const auto &pwhs = geo_divs[g].polygons_with_holes();
    for (std::size_t p = 0; p < pwhs.size(); ++p) {
      rings.emplace_back(
        &pwhs_dens[p].outer_boundary(),
        &pwhs[p].outer_boundary());
      auto h_dens = pwhs_dens[p].holes().begin();
      for (const auto &h : pwhs[p].holes()) {
        rings.emplace_back(&*h_dens, &h);
        ++h_dens;
      }
    }
  }
  assert(rings.size() == inserted.size());

  std::vector<bool> fails(rings.size(), false);
  std::vector<bool> changed(rings.size(), false);
  std::vector<std::size_t> to_check(rings.size());
  std::iota(to_check.begin(), to_check.end(), std::size_t{0});
  std::unordered_set<Point> removed;
  while (!to_check.empty()) {
    std::unordered_set<Point> drop;
    for (const std::size_t r : to_check) {
      const Polygon &ring = *rings[r].first;
      const Polygon &before = *rings[r].second;
      fails[r] = !ring.is_simple();
      const bool repair =
        (fails[r] && before.is_simple()) ||
        (changed[r] && !is_simple_with_exact_predicates(ring) &&
         is_simple_with_exact_predicates(before));
      if (!repair)
        continue;
      std::unordered_set<Point> found;
      densification_points_on_intersecting_segments(ring, inserted[r], found);
      if (found.empty()) {
        for (std::size_t k = 0; k < ring.size(); ++k) {
          if (inserted[r][k])
            found.insert(ring[k]);
        }
      }
      drop.insert(found.begin(), found.end());
    }

    if (drop.empty())
      break;

    to_check.clear();
    for (std::size_t r = 0; r < rings.size(); ++r) {
      Polygon &ring = *rings[r].first;
      std::vector<bool> &ins = inserted[r];
      Polygon kept;
      std::vector<bool> kept_ins;
      for (std::size_t k = 0; k < ring.size(); ++k) {
        if (ins[k] && drop.count(ring[k]) > 0) {
          removed.insert(ring[k]);
        } else {
          kept.push_back(ring[k]);
          kept_ins.push_back(ins[k]);
        }
      }
      if (kept.size() != ring.size()) {
        ring = std::move(kept);
        ins = std::move(kept_ins);
        changed[r] = true;
        to_check.push_back(r);
      }
    }
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

void InsetState::densify_geo_divs_using_delaunay_t()
{
  timer.start("Densification");

  std::cerr << "Densifying" << std::endl;
  size_t n_points_before = n_points();
  std::cerr << "Num points before densification: " << n_points_before
            << std::endl;

  const auto &triangles = triang_.triangles();

  std::vector<std::array<uint32_t, 4>> edges_raw;
  edges_raw.reserve(triangles.size() * 3);

  auto push = [&](const auto &p, const auto &q) {
    if (q.x() < p.x() || (q.x() == p.x() && q.y() < p.y())) {
      edges_raw.push_back({q.x(), q.y(), p.x(), p.y()});
    } else {
      edges_raw.push_back({p.x(), p.y(), q.x(), q.y()});
    }
  };

  for (const auto &tri : triangles) {
    push(tri.vertices[0], tri.vertices[1]);
    push(tri.vertices[1], tri.vertices[2]);
    push(tri.vertices[2], tri.vertices[0]);
  }

  std::sort(edges_raw.begin(), edges_raw.end());
  edges_raw.erase(
    std::unique(edges_raw.begin(), edges_raw.end()),
    edges_raw.end());

  // Pre-classify edges (H/V/OTHER), so that crossings with horizontal and
  // vertical edges can be placed exactly on their lattice lines
  std::vector<EdgeKind> edge_kind;
  edge_kind.reserve(edges_raw.size());
  for (const auto &s : edges_raw)
    edge_kind.push_back(classify_edge(s[0], s[1], s[2], s[3]));

  std::vector<Segment> edge_vec;
  edge_vec.reserve(edges_raw.size());
  for (auto const &e : edges_raw) {
    edge_vec.emplace_back(Point(e[0], e[1]), Point(e[2], e[3]));
  }

  AABB_Tree tree(edge_vec.begin(), edge_vec.end());
  tree.build();

  // We have lots of shared boundaries, so we will cache the results of
  // densification points for each edge, to avoid recomputing them twice
  SegCache cache;
  cache.reserve(n_points_before);

  std::vector<GeoDiv> geodivs_dens;
  geodivs_dens.reserve(geo_divs_.size());

  // For each ring, in the order of the rings in geodivs_dens, whether
  // densification inserted each point
  std::vector<std::vector<bool>> inserted;

  for (const auto &gd : geo_divs_) {
    GeoDiv gd_dens(gd.id());

    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto &outer = pwh.outer_boundary();
      Polygon outer_dens;
      std::vector<bool> inserted_in_outer;
      outer_dens.reserve(outer.size() * 5);

      const std::size_t outer_sz = outer.size();

      for (std::size_t i = 0; i < outer_sz; ++i) {
        const Point a = outer[i];
        const Point b = (i + 1 == outer_sz) ? outer[0] : outer[i + 1];

        const std::vector<Point> seg_pts = densification_points_with_edge_tree(
          a,
          b,
          tree,
          edge_vec,
          edge_kind,
          &cache);

        append_densified_edge(outer_dens, inserted_in_outer, seg_pts);
      }
      inserted.push_back(std::move(inserted_in_outer));

      std::vector<Polygon> holes_v_dens;
      holes_v_dens.reserve(pwh.number_of_holes());

      for (auto const &h : pwh.holes()) {
        Polygon hole_dens;
        std::vector<bool> inserted_in_hole;
        hole_dens.reserve(h.size() * 2);

        const std::size_t h_sz = h.size();
        for (std::size_t j = 0; j < h_sz; ++j) {
          const Point c = h[j];
          const Point d = (j + 1 == h_sz) ? h[0] : h[j + 1];

          const std::vector<Point> seg_pts =
            densification_points_with_edge_tree(
              c,
              d,
              tree,
              edge_vec,
              edge_kind,
              &cache);

          append_densified_edge(hole_dens, inserted_in_hole, seg_pts);
        }

        inserted.push_back(std::move(inserted_in_hole));
        holes_v_dens.emplace_back(std::move(hole_dens));
      }

      const Polygon_with_holes pwh_dens(
        std::move(outer_dens),
        holes_v_dens.begin(),
        holes_v_dens.end());
      gd_dens.push_back(std::move(pwh_dens));
    }

    geodivs_dens.emplace_back(std::move(gd_dens));
  }

  const bool all_simple = drop_densification_points_that_break_simplicity(
    geodivs_dens,
    geo_divs_,
    inserted);
  geo_divs_ = std::move(geodivs_dens);

  std::cerr << "Num points after densification: " << n_points() << std::endl;

  // Every ring has passed the test of is_simple(), unless all_simple is
  // false; is_simple() then reports the first ring that fails
  if (!all_simple)
    is_simple(__func__);
  timer.stop("Densification");
}
