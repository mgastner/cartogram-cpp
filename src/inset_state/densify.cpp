#include "inset_state.hpp"
#include "round_point.hpp"
#include <CGAL/AABB_segment_primitive_2.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/intersections.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <unordered_map>

using AABB_Prim =
  CGAL::AABB_segment_primitive_2<Scd, std::vector<Segment>::const_iterator>;
using AABB_Traits = CGAL::AABB_traits_2<Scd, AABB_Prim>;
using AABB_Tree = CGAL::AABB_tree<AABB_Traits>;

using SegCache = std::unordered_map<Segment, std::vector<Point>>;

// Edge classification of AABB tree segments
// AABB tree primitives give the idex of the segment that of the vector we
// passed to create the tree, so we can classify edges based on their index in
// that vector, and use fast-paths for common cases, avoiding expensive CGAL
// intersection calls; This is particularly useful for our case because most of
// our edges are axis-aligned or 45° diagonal, and we can handle them with
// simple arithmetic.

// H    (0): horizontal edge, y1 == y2
// V    (1): vertical edge,   x1 == x2
// DPOS (2): 45° diagonal with positive slope (|dx|==|dy| and dx*dy > 0)  ↗ / ↙
// DNEG (3): 45° diagonal with negative slope (|dx|==|dy| and dx*dy < 0)  ↖ / ↘
// OTHER(4): all other edges (non-axis-aligned, non-45°, or degenerate)
enum class EdgeKind : uint8_t { H = 0, V = 1, DPOS = 2, DNEG = 3, OTHER = 4 };

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

  const int64_t dx = int64_t(x1) - int64_t(x0);
  const int64_t dy = int64_t(y1) - int64_t(y0);

  if (dy == dx)
    return EdgeKind::DPOS;  // slope +1
  if (dy == -dx)
    return EdgeKind::DNEG;  // slope -1

  return EdgeKind::OTHER;
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

  // Param t along [pt1, pt2] for ordering hits
  const double ax = pt1.x(), ay = pt1.y();
  const double bx = pt2.x(), by = pt2.y();
  const double dax = bx - ax, day = by - ay;
  const double L2 = dax * dax + day * day;

  auto t_of = [&](const Point &p) -> double {
    if (L2 == 0.0)
      return 0.0;
    return ((p.x() - ax) * dax + (p.y() - ay) * day) / L2;
  };

  // Collect primitives that intersect the query (not just bbox hits)
  std::vector<AABB_Tree::Primitive_id> prims;
  prims.reserve(8);
  tree.all_intersected_primitives(query, std::back_inserter(prims));

  // Always include endpoints
  std::vector<std::pair<double, Point>> hits;
  hits.reserve(prims.size() * 2 + 2);
  hits.emplace_back(0.0, pt1);
  hits.emplace_back(1.0, pt2);

  for (const auto &pid : prims) {

    // Primitive_id is an iterator into edge_vec
    const std::ptrdiff_t idx = pid - edge_vec.begin();

    assert(idx >= 0 && static_cast<size_t>(idx) < edge_vec.size());

    const Segment &E = *pid;
    const EdgeKind kind = edge_kind[static_cast<size_t>(idx)];

    const double ex0 = E.source().x(), ey0 = E.source().y();
    const double ex1 = E.target().x(), ey1 = E.target().y();

    auto in_closed01 = [&](double t) {
      return less_than_equal(0.0, t) && less_than_equal(t, 1.0);
    };

    auto within_closed = [&](double v, double a, double b) {
      const double lo = std::min(a, b);
      const double hi = std::max(a, b);
      return less_than_equal(lo, v) && less_than_equal(v, hi);
    };

    // Here we use the edge classification to handle the most common cases
    // with simple arithmetic, avoiding expensive CGAL intersection calls
    switch (kind) {

    case EdgeKind::H: {
      const double y0 = ey0;  // == ey1
      const double dy = day;  // query dy
      if (almost_equal(dy, 0.0)) {
        // Query horizontal: either collinear or disjoint
        if (almost_equal(ay, y0)) {
          // Collinear overlap: push both endpoints of edge (clipped to [0,1])
          const Point p0(ex0, y0), p1(ex1, y0);
          const double t0 = t_of(p0), t1 = t_of(p1);
          if (in_closed01(t0))
            hits.emplace_back(t0, p0);
          if (in_closed01(t1))
            hits.emplace_back(t1, p1);
        }
        // else disjoint: nothing to add
      } else {
        const double t = (y0 - ay) / dy;
        if (!in_closed01(t)) {
          // intersection with infinite line lies outside the segment
        } else {
          const double x = ax + t * dax;
          if (within_closed(x, ex0, ex1))
            hits.emplace_back(std::clamp(t, 0.0, 1.0), Point(x, y0));
        }
      }
      break;
    }

    case EdgeKind::V: {
      const double x0 = ex0;  // == ex1
      const double dx = dax;  // query dx
      if (almost_equal(dx, 0.0)) {
        // Query vertical: collinear or disjoint
        if (almost_equal(ax, x0)) {
          const Point p0(x0, ey0), p1(x0, ey1);
          const double t0 = t_of(p0), t1 = t_of(p1);
          if (in_closed01(t0))
            hits.emplace_back(t0, p0);
          if (in_closed01(t1))
            hits.emplace_back(t1, p1);
        }
      } else {
        const double t = (x0 - ax) / dx;
        if (!in_closed01(t)) {
          // outside segment
        } else {
          const double y = ay + t * day;
          if (within_closed(y, ey0, ey1))
            hits.emplace_back(std::clamp(t, 0.0, 1.0), Point(x0, y));
        }
      }
      break;
    }

    case EdgeKind::DPOS: {
      // Edge: y - ey0 = +(x - ex0)  =>  y - x = c
      const double c = ey0 - ex0;
      const double denom = (day - dax);  // (qy - qx)
      if (almost_equal(denom, 0.0)) {
        // Query also slope +1; collinear iff same c
        const double c_q = (ay - ax);
        if (almost_equal(c_q, c)) {
          const Point p0(ex0, ey0), p1(ex1, ey1);
          const double t0 = t_of(p0), t1 = t_of(p1);
          if (in_closed01(t0))
            hits.emplace_back(t0, p0);
          if (in_closed01(t1))
            hits.emplace_back(t1, p1);
        }
      } else {
        const double t = (c - (ay - ax)) / denom;
        if (in_closed01(t)) {
          const double x = ax + t * dax;
          const double y = ay + t * day;
          if (within_closed(x, ex0, ex1) && within_closed(y, ey0, ey1))
            hits.emplace_back(std::clamp(t, 0.0, 1.0), Point(x, y));
        }
      }
      break;
    }

    case EdgeKind::DNEG: {
      // Edge: y - ey0 = -(x - ex0)  =>  y + x = c
      const double c = ey0 + ex0;
      const double denom = (day + dax);  // (qy + qx)
      if (almost_equal(denom, 0.0)) {
        const double c_q = (ay + ax);
        if (almost_equal(c_q, c)) {
          const Point p0(ex0, ey0), p1(ex1, ey1);
          const double t0 = t_of(p0), t1 = t_of(p1);
          if (in_closed01(t0))
            hits.emplace_back(t0, p0);
          if (in_closed01(t1))
            hits.emplace_back(t1, p1);
        }
      } else {
        const double t = (c - (ay + ax)) / denom;
        if (in_closed01(t)) {
          const double x = ax + t * dax;
          const double y = ay + t * day;
          if (within_closed(x, ex0, ex1) && within_closed(y, ey0, ey1))
            hits.emplace_back(std::clamp(t, 0.0, 1.0), Point(x, y));
        }
      }
      break;
    }

    case EdgeKind::OTHER: {
      CGAL::Object obj = CGAL::intersection(query, E);
      Point ip;
      if (CGAL::assign(ip, obj)) {
        const double t = t_of(ip);
        if (in_closed01(t))
          hits.emplace_back(t, ip);
      } else {
        Segment seg;
        if (CGAL::assign(seg, obj)) {
          const Point p0 = seg.source();
          const Point p1 = seg.target();
          const double t0 = t_of(p0);
          const double t1 = t_of(p1);
          if (in_closed01(t0))
            hits.emplace_back(t0, p0);
          if (in_closed01(t1))
            hits.emplace_back(t1, p1);
        }
      }
      break;
    }
    default:
      // Should not happen
      assert(false);
      break;
    }
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

  // Very important step: Pre-classify edges (H/V/diag/OTHER)
  // This allows us to later handle intersection computation of the most common
  // cases with simple arithmetic
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

  for (const auto &gd : geo_divs_) {
    GeoDiv gd_dens(gd.id());

    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto &outer = pwh.outer_boundary();
      Polygon outer_dens;
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

        if (seg_pts.size() > 1)
          outer_dens.insert(
            outer_dens.end(),
            seg_pts.begin(),
            seg_pts.end() - 1);
      }

      std::vector<Polygon> holes_v_dens;
      holes_v_dens.reserve(pwh.number_of_holes());

      for (auto const &h : pwh.holes()) {
        Polygon hole_dens;
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

          if (seg_pts.size() > 1)
            hole_dens.insert(
              hole_dens.end(),
              seg_pts.begin(),
              seg_pts.end() - 1);
        }

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

  geo_divs_ = std::move(geodivs_dens);

  std::cerr << "Num points after densification: " << n_points() << std::endl;
  is_simple(__func__);
  timer.stop("Densification");
}
