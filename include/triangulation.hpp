#pragma once
#include "cgal_typedef.hpp"
#include "quadtree_corner.hpp"
#include "round_point.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/enum.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPoint = EPICK::Point_2;

template <class QuadtreeLocator, class Projection> class Triangulation
{
public:
  struct Triangle {
    std::array<QuadtreeCorner, 3> vertices;  // grid corners (x, y)
  };

  using TriangleIt = typename std::vector<Triangle>::iterator;

  // We first build per Quadtree leaf triangulation (we choose the
  // triangulation that maximizes the minimum angle in the projected space).
  // However, if the triangulation triangle flips in the projected space that
  // implies our blur width is too small. We return false in that case
  bool build(const QuadtreeLocator *qt_locator, const Projection *proj_data)
  {
    clear();

    // Update pointers to the current data
    qt_locator_ = qt_locator;
    proj_data_ = proj_data;

    const auto &leaves = qt_locator_->leaves();
    triangles_.reserve(leaves.size() * 4);
    leaf_to_triangle_idxs_.resize(leaves.size() * 4);

    // The leaves from the qt_locator_ only contains the bottom-left corner
    // (x, y) and the size of the leaf. In order to triangulate the leaf, we
    // also need the potential midpoints of the edges (if the quadtree cells
    // besides the leaf have higher depth)

    // Now for each leaf, we first construct the polygon of the leaf including
    // the potential midpoints of the edges. Then we triangulate the polygon
    // to maximize the minimum angle in the projected space
    for (const auto &leaf : leaves) {
      const uint32_t x = leaf.x;
      const uint32_t y = leaf.y;
      const uint32_t s = leaf.size;

      // Collect boundary in CCW: bl -> br -> tr -> tl (+ midpoints as needed)
      std::vector<QuadtreeCorner> leaf_pts;
      leaf_pts.reserve(8);

      // bottom edge (bl -> br)
      leaf_pts.push_back({x, y});

      if (s > 1) {
        const uint32_t x_mid = x + (s >> 1);
        const uint32_t y_mid = y;
        if (proj_data_->is_valid_corner(x_mid, y_mid)) {
          leaf_pts.push_back({x_mid, y_mid});  // midpoint
        }
      }
      leaf_pts.push_back({x + s, y});

      // right edge (br -> tr)
      if (s > 1) {
        const uint32_t x_mid = x + s;
        const uint32_t y_mid = y + (s >> 1);
        if (proj_data_->is_valid_corner(x_mid, y_mid)) {
          leaf_pts.push_back({x_mid, y_mid});  // midpoint
        }
      }
      leaf_pts.push_back({x + s, y + s});

      // top edge (tr -> tl)
      if (s > 1) {
        const uint32_t x_mid = x + (s >> 1);
        const uint32_t y_mid = y + s;
        if (proj_data_->is_valid_corner(x_mid, y_mid)) {
          leaf_pts.push_back({x_mid, y_mid});  // midpoint
        }
      }
      leaf_pts.push_back({x, y + s});

      // left edge (tl -> bl)
      if (s > 1) {
        const uint32_t x_mid = x;
        const uint32_t y_mid = y + (s >> 1);
        if (proj_data_->is_valid_corner(x_mid, y_mid)) {
          leaf_pts.push_back({x_mid, y_mid});  // midpoint
        }
      }

      std::vector<Point> leaf_pts_proj;
      leaf_pts_proj.reserve(leaf_pts.size());
      for (size_t i = 0; i < leaf_pts.size(); ++i) {
        const auto &pt = leaf_pts[i];
        leaf_pts_proj.push_back(proj_data_->get(pt.x(), pt.y()));
      }

      // Triangulate polygon to maximize minimum angle in projected space
      // We return false from the triangulation in case any triangle of the
      // optimal triangulation of the leaf flips in the projected space
      if (!triangulate_max_min_angle(leaf_pts, leaf_pts_proj, x, y)) {
        return false;
      }
    }

    return true;
  }

  // Internally, we use the CGAL EPICK-based containment to locate the triangle
  // that contains the point p. The benefit of this is that we do not have to
  // worry about choosing the wrong triangle because of some bad choice of
  // EPSILON With EPICK, all predicates such as orientation, contains, etc. are
  // exact, so we can use them directly without worrying about numerical issues
  TriangleIt locate(const Point &p)
  {
    auto contains = [&](const Triangle &T) -> bool {
      const EPoint A(T.vertices[0].x(), T.vertices[0].y());
      const EPoint B(T.vertices[1].x(), T.vertices[1].y());
      const EPoint C(T.vertices[2].x(), T.vertices[2].y());
      const EPoint P(p.x(), p.y());

      const auto orientation = CGAL::orientation(A, B, C);

      assert(
        orientation != CGAL::COLLINEAR &&
        "Collinear triangles found in triangulation.");

      if (orientation == CGAL::LEFT_TURN) {
        return CGAL::orientation(A, B, P) != CGAL::RIGHT_TURN &&
               CGAL::orientation(B, C, P) != CGAL::RIGHT_TURN &&
               CGAL::orientation(C, A, P) != CGAL::RIGHT_TURN;
      } else {  // RIGHT_TURN
        return CGAL::orientation(A, B, P) != CGAL::LEFT_TURN &&
               CGAL::orientation(B, C, P) != CGAL::LEFT_TURN &&
               CGAL::orientation(C, A, P) != CGAL::LEFT_TURN;
      }
    };

    // Fast path -> check if the last located triangle contains the point
    // At the later stages, it is highly likely that many points will be within
    // the same triangle, so we can avoid the full `locate` cost by just doing
    // a low cost containment check in the last located triangle
    if (last_locate_triangle_idx_ != UINT32_MAX) {
      if (contains(triangles_[last_locate_triangle_idx_]))
        return triangles_.begin() + last_locate_triangle_idx_;
    }

    // First find the leaf from quadtree that contains the point p
    // Then check all triangles in that leaf (should be quite fast because each
    // leaf only has a handful of triangles)
    const auto leaf = qt_locator_->locate(p);
    const uint32_t key = proj_data_->offset(leaf.x, leaf.y);
    assert(key < leaf_to_triangle_idxs_.size());

    const auto &triangle_indices = leaf_to_triangle_idxs_[key];
    assert(!triangle_indices.empty());

    for (const uint32_t idx : triangle_indices) {
      if (contains(triangles_[idx])) {
        last_locate_triangle_idx_ = idx;
        return triangles_.begin() + idx;
      }
    }

    assert(false && "No triangle found for the point in the quadtree leaf");
    return triangles_.end();
  }

  [[nodiscard]] const std::vector<Triangle> &triangles() const
  {
    return triangles_;
  }

private:
  const QuadtreeLocator *qt_locator_{};
  const Projection *proj_data_{};

  std::vector<Triangle> triangles_;
  std::vector<std::vector<uint32_t>> leaf_to_triangle_idxs_;

  mutable uint32_t last_locate_triangle_idx_{UINT32_MAX};

  void clear()
  {
    triangles_.clear();
    leaf_to_triangle_idxs_.clear();
    last_locate_triangle_idx_ = UINT32_MAX;
  }

  // `leaf_pts` are the original-space vertices (integral values given by
  // x(), y()) and `leaf_pts_proj` are the projected-space vertices (floating
  // value given by x(), y()). The `leaf_x` and `leaf_y`
  // are the coordinates of the leaf in the quadtree. Later we store the
  // triangles in the `triangles_` vector and the indices of the triangles in
  // the `leaf_to_triangle_idxs_` vector, where the key is the offset of the
  // leaf in the projection data (which can be accessed via
  // `proj_data_->offset(leaf_x, leaf_y)`)
  bool triangulate_max_min_angle(
    const std::vector<QuadtreeCorner> &leaf_pts,
    const std::vector<Point> &leaf_pts_proj,
    uint32_t leaf_x,
    uint32_t leaf_y)
  {
    const uint32_t n = static_cast<uint32_t>(leaf_pts.size());

    // Given that we are triangulating a square with possibly 0 to 4 midpoints
    // in addtion to the 4 corners, the following condition must be true
    assert(n >= 4 && n <= 8);

    // Sanity check
    assert(leaf_pts_proj.size() == n);

    auto angle_at =
      [&](const Point &p, const Point &q, const Point &r) -> double {
      const double ux = p.x() - q.x(), uy = p.y() - q.y();
      const double vx = r.x() - q.x(), vy = r.y() - q.y();
      return std::atan2(
        std::abs(ux * vy - uy * vx),
        ux * vx + uy * vy);  // in [0, π]
    };

    auto compute_min_angle =
      [&](const Point &a, const Point &b, const Point &c) {
        const double a0 = angle_at(b, a, c);
        const double a1 = angle_at(a, b, c);
        const double a2 = angle_at(a, c, b);
        return std::min({a0, a1, a2});
      };

    // Since the coordinates are all integral, we can do the check directly
    // without worrying about floating point precision issues
    auto axis_colinear = [](
                           const QuadtreeCorner &a,
                           const QuadtreeCorner &b,
                           const QuadtreeCorner &c) -> bool {
      const auto ax = a.x(), ay = a.y();
      const auto bx = b.x(), by = b.y();
      const auto cx = c.x(), cy = c.y();
      return (ax == bx && bx == cx) || (ay == by && by == cy);
    };

    // Here we use dynamic programming to find the triangulation that maximizes
    // the minimum angle in the projected space

    // DP: dp[i][j] = best achievable minimal angle (radians) on chain i..j
    // (i<j)
    std::vector<std::vector<double>> dp(n, std::vector<double>(n, 0.0));

    // For each final value of chain i..j in dp[i][j], cut[i][j] where should
    // the chain be cut at k, so that min(dp[i][k], min_angle(i, k, j),
    // dp[k][j]) = dp[i][j]
    // This is useful to construct the triangulation later
    std::vector<std::vector<int>> cut(n, std::vector<int>(n, -1));

    // Base case with triangle with two points
    for (size_t i = 0; i + 1 < n; ++i)
      dp[i][i + 1] = std::numeric_limits<double>::infinity();

    for (size_t intv = 2; intv < n; ++intv) {
      for (size_t i = 0; i + intv < n; ++i) {
        size_t j = i + intv;

        // Now need to find the best "cut" k that will maximize dp[i][j]
        // Notice that values dp[i][k], and dp[k][j] are already computed
        // previously
        double best = 0.0;
        int bestk = -1;
        for (size_t k = i + 1; k <= j - 1; ++k) {

          // We do not want to consider triangles that are not triangles at all
          if (axis_colinear(leaf_pts[i], leaf_pts[k], leaf_pts[j]))
            continue;

          const double min_proj_angle = compute_min_angle(
            leaf_pts_proj[i],
            leaf_pts_proj[k],
            leaf_pts_proj[j]);

          const double cand = std::min({dp[i][k], dp[k][j], min_proj_angle});
          if (less_than(best, cand)) {
            best = cand;
            bestk = static_cast<int>(k);
          }
        }

        dp[i][j] = best;
        cut[i][j] = bestk;
      }
    }

    // Reconstruct triangles for the optimal solution at (dp[0][n-1])

    bool okay =
      true;  // If any triangle flips in the projected space, we return false
    size_t leaf_offset = proj_data_->offset(leaf_x, leaf_y);
    auto construct = [&](auto &&self, size_t i, size_t j) -> void {
      if (j <= i + 1)
        return;
      assert(cut[i][j] >= 0);
      size_t k = static_cast<size_t>(cut[i][j]);
      assert(k >= 0);

      CGAL::Orientation orien_ori = CGAL::orientation(
        Point(leaf_pts[i]),
        Point(leaf_pts[k]),
        Point(leaf_pts[j]));
      CGAL::Orientation orien_proj = CGAL::orientation(
        leaf_pts_proj[i],
        leaf_pts_proj[k],
        leaf_pts_proj[j]);

      // Check if triangles flips or either of them is collinear
      if (orien_proj == CGAL::COLLINEAR || orien_ori != orien_proj) {
        okay = false;
      }

      Triangle T{{leaf_pts[i], leaf_pts[k], leaf_pts[j]}};
      const uint32_t triangle_idx = static_cast<uint32_t>(triangles_.size());
      triangles_.push_back(T);
      leaf_to_triangle_idxs_[leaf_offset].push_back(triangle_idx);

      self(self, i, static_cast<uint32_t>(k));
      self(self, static_cast<uint32_t>(k), j);
    };

    construct(construct, 0, n - 1);
    return okay;
  }
};
