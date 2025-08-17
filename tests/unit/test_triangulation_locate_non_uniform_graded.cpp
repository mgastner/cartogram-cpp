#define BOOST_TEST_MODULE test_triangulation_locate_non_uniform_graded
#include "cgal_typedef.hpp"
#include "quadtree_corner.hpp"
#include "round_point.hpp"
#include "triangulation.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <algorithm>
#include <array>
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPoint = EPICK::Point_2;

static constexpr uint32_t GRID_W = 8;
static constexpr uint32_t GRID_H = 8;

struct NonUniformQuadtreeLocator {
  struct Leaf {
    uint32_t x, y, size;
  };

  NonUniformQuadtreeLocator() = default;

  explicit NonUniformQuadtreeLocator(
    std::vector<Leaf> leaves,
    uint32_t w = GRID_W,
    uint32_t h = GRID_H)
      : w_(w), h_(h), leaves_(std::move(leaves))
  {
    build_maps_();
  }

  const std::vector<Leaf> &leaves() const
  {
    return leaves_;
  }

  Leaf locate(const Point &p) const
  {
    const double xd = CGAL::to_double(p.x());
    const double yd = CGAL::to_double(p.y());
    const int ix = std::max(
      0,
      std::min<int>(
        static_cast<int>(std::floor(xd)),
        static_cast<int>(w_) - 1));
    const int iy = std::max(
      0,
      std::min<int>(
        static_cast<int>(std::floor(yd)),
        static_cast<int>(h_) - 1));
    const int idx = cell_leaf_idx_[ix][iy];
    BOOST_REQUIRE(idx >= 0);
    return leaves_[static_cast<size_t>(idx)];
  }

  uint32_t width() const
  {
    return w_;
  }
  uint32_t height() const
  {
    return h_;
  }

  size_t count_midpoints(const Leaf &L) const
  {
    return edge_mask_smaller(L).count();
  }

  // Bitset over 4 edges (0: bottom, 1: right, 2: top, 3: left) where neighbor
  // is smaller (s/2).
  std::bitset<4> edge_mask_smaller(const Leaf &L) const
  {
    std::bitset<4> mask;
    const uint32_t x = L.x, y = L.y, s = L.size;
    if (s > 1) {
      if (y > 0) {
        mask.set(0, size_at_(x, y - 1) < s);
      }
      if (x + s < w_) {
        mask.set(1, size_at_(x + s, y) < s);
      }
      if (y + s < h_) {
        mask.set(2, size_at_(x, y + s) < s);
      }
      if (x > 0) {
        mask.set(3, size_at_(x - 1, y) < s);
      }
    }
    return mask;
  }

  // Layout A: mixture of 4x4, 2x2, 1x1 (balanced).
  // - TR: one 4x4 at (4,4)
  // - TL: four 2x2 over [0,4)x[4,8)
  // - BR: four 2x2 over [4,8)x[0,4)
  // - BL: mostly 2x2 with a 2x2 block of 1x1 at (0..2,0..2)
  static NonUniformQuadtreeLocator make_layout_A()
  {
    std::vector<Leaf> L;

    // TR 4x4
    L.push_back({4, 4, 4});

    // TL 2x2
    L.push_back({0, 4, 2});
    L.push_back({2, 4, 2});
    L.push_back({0, 6, 2});
    L.push_back({2, 6, 2});

    // BR 2x2
    L.push_back({4, 0, 2});
    L.push_back({6, 0, 2});
    L.push_back({4, 2, 2});
    L.push_back({6, 2, 2});

    // BL 1x1 patch + surrounding 2x2
    L.push_back({0, 0, 1});
    L.push_back({1, 0, 1});
    L.push_back({0, 1, 1});
    L.push_back({1, 1, 1});
    L.push_back({2, 0, 2});
    L.push_back({0, 2, 2});
    L.push_back({2, 2, 2});

    return NonUniformQuadtreeLocator(std::move(L));
  }

  // Layout B: all 2x2, except a central 2x2 split to four 1x1 (balanced).
  static NonUniformQuadtreeLocator make_layout_B()
  {
    std::vector<Leaf> L;

    // 2x2 tiling
    for (uint32_t y = 0; y < GRID_H; y += 2)
      for (uint32_t x = 0; x < GRID_W; x += 2)
        L.push_back({x, y, 2});

    // Refine the 2x2 leaf at (2,2) to 1x1
    // Remove the {2,2,2} entry and add 4x{1,1}
    L.erase(
      std::remove_if(
        L.begin(),
        L.end(),
        [](const Leaf &a) {
          return a.x == 2 && a.y == 2 && a.size == 2;
        }),
      L.end());
    L.push_back({2, 2, 1});
    L.push_back({3, 2, 1});
    L.push_back({2, 3, 1});
    L.push_back({3, 3, 1});
    return NonUniformQuadtreeLocator(std::move(L));
  }

private:
  uint32_t w_{GRID_W}, h_{GRID_H};
  std::vector<Leaf> leaves_;

  int cell_leaf_idx_[GRID_W][GRID_H]{};  // index into leaves_
  uint32_t cell_leaf_size_[GRID_W][GRID_H]{};  // size of the covering leaf

  void build_maps_()
  {
    for (uint32_t y = 0; y < h_; ++y)
      for (uint32_t x = 0; x < w_; ++x)
        cell_leaf_idx_[x][y] = -1, cell_leaf_size_[x][y] = 0;

    // Paint leaf index + size over unit cells
    for (size_t i = 0; i < leaves_.size(); ++i) {
      const auto &L = leaves_[i];
      for (uint32_t yy = L.y; yy < L.y + L.size; ++yy)
        for (uint32_t xx = L.x; xx < L.x + L.size; ++xx) {
          BOOST_REQUIRE(cell_leaf_idx_[xx][yy] == -1);  // no overlaps
          cell_leaf_idx_[xx][yy] = static_cast<int>(i);
          cell_leaf_size_[xx][yy] = L.size;
        }
    }

    // Full coverage
    for (uint32_t y = 0; y < h_; ++y)
      for (uint32_t x = 0; x < w_; ++x)
        BOOST_REQUIRE(cell_leaf_idx_[x][y] >= 0);
  }

  uint32_t size_at_(uint32_t cell_x, uint32_t cell_y) const
  {
    BOOST_REQUIRE(cell_x < w_ && cell_y < h_);
    return cell_leaf_size_[cell_x][cell_y];
  }
};

static inline uint64_t pack_xy(uint32_t x, uint32_t y)
{
  return (uint64_t(x) << 32) | uint64_t(y);
}

template <class Deform> struct BalancedProjection {
  BalancedProjection(const NonUniformQuadtreeLocator *qt, Deform deform)
      : qt_(qt), deform_(std::move(deform)), w_(qt->width()), h_(qt->height())
  {
    const auto &leaves = qt_->leaves();
    for (size_t i = 0; i < leaves.size(); ++i) {
      const auto &L = leaves[i];
      offset_map_[pack_xy(L.x, L.y)] = static_cast<uint32_t>(i * 4);
    }

    for (const auto &L : leaves) {
      const uint32_t x = L.x, y = L.y, s = L.size;
      if (s <= 1)
        continue;
      const auto mask = qt_->edge_mask_smaller(L);
      if (mask.test(0))
        midpoints_.insert(pack_xy(x + (s >> 1), y));
      if (mask.test(1))
        midpoints_.insert(pack_xy(x + s, y + (s >> 1)));
      if (mask.test(2))
        midpoints_.insert(pack_xy(x + (s >> 1), y + s));
      if (mask.test(3))
        midpoints_.insert(pack_xy(x, y + (s >> 1)));
    }
  }

  bool is_valid_corner(uint32_t x, uint32_t y) const
  {
    return midpoints_.find(pack_xy(x, y)) != midpoints_.end();
  }

  uint32_t offset(uint32_t x, uint32_t y) const
  {
    auto it = offset_map_.find(pack_xy(x, y));
    BOOST_REQUIRE(it != offset_map_.end());
    return it->second;
  }

  Point get(uint32_t x, uint32_t y) const
  {
    const auto p = deform_(
      static_cast<double>(x),
      static_cast<double>(y),
      static_cast<double>(w_),
      static_cast<double>(h_));
    return Point(p.first, p.second);
  }

  size_t num_unique_corners() const
  {
    return (w_ + 1) * (h_ + 1);
  }

private:
  const NonUniformQuadtreeLocator *qt_{};
  Deform deform_;
  uint32_t w_, h_;
  std::unordered_map<uint64_t, uint32_t> offset_map_;
  std::unordered_set<uint64_t> midpoints_;
};

struct ShearX {
  double k{};
  std::pair<double, double> operator()(
    double x,
    double y,
    double /*w*/,
    double /*h*/) const
  {
    return {x + k * y, y};
  }
};

struct Rotate {
  double theta{};
  std::pair<double, double> operator()(double x, double y, double w, double h)
    const
  {
    const double cx = 0.5 * w, cy = 0.5 * h;
    const double c = std::cos(theta), s = std::sin(theta);
    const double X = x - cx, Y = y - cy;
    return {c * X - s * Y + cx, s * X + c * Y + cy};
  }
};

struct RadialInward {
  double alpha{};  // in (0,1)
  std::pair<double, double> operator()(double x, double y, double w, double h)
    const
  {
    const double cx = 0.5 * w, cy = 0.5 * h;
    const double s = (1.0 - alpha);
    return {s * (x - cx) + cx, s * (y - cy) + cy};
  }
};

struct Wavy {
  double eps{};
  std::pair<double, double> operator()(
    double x,
    double y,
    double /*w*/,
    double /*h*/) const
  {
    return {x + eps * std::sin(0.7 * y), y + eps * std::sin(0.5 * x)};
  }
};

struct MirrorX {
  std::pair<double, double> operator()(
    double x,
    double /*y*/,
    double w,
    double h) const
  {
    (void)h;
    return {
      w - x,
      /*y unchanged*/ 0.0};
  }
};

struct MirrorXAdapter {
  std::pair<double, double> operator()(double x, double y, double w, double h)
    const
  {
    (void)h;
    return {w - x, y};
  }
};

template <class TriangulationT>
static bool triangle_contains(
  const typename TriangulationT::Triangle &T,
  const Point &p)
{
  const EPoint A(T.vertices[0].x(), T.vertices[0].y());
  const EPoint B(T.vertices[1].x(), T.vertices[1].y());
  const EPoint C(T.vertices[2].x(), T.vertices[2].y());
  const EPoint P(p.x(), p.y());

  const auto orientation = CGAL::orientation(A, B, C);
  BOOST_REQUIRE(orientation != CGAL::COLLINEAR);

  if (orientation == CGAL::LEFT_TURN) {
    return CGAL::orientation(A, B, P) != CGAL::RIGHT_TURN &&
           CGAL::orientation(B, C, P) != CGAL::RIGHT_TURN &&
           CGAL::orientation(C, A, P) != CGAL::RIGHT_TURN;
  } else {
    return CGAL::orientation(A, B, P) != CGAL::LEFT_TURN &&
           CGAL::orientation(B, C, P) != CGAL::LEFT_TURN &&
           CGAL::orientation(C, A, P) != CGAL::LEFT_TURN;
  }
}

// Return allowed boundary vertices for a leaf: 4 corners + present hanging
// nodes
static std::vector<std::pair<uint32_t, uint32_t>> leaf_allowed_vertices(
  const NonUniformQuadtreeLocator &qt,
  const NonUniformQuadtreeLocator::Leaf &L)
{
  std::vector<std::pair<uint32_t, uint32_t>> V;
  const uint32_t x = L.x, y = L.y, s = L.size;
  V.reserve(8);
  V.emplace_back(x, y);
  V.emplace_back(x + s, y);
  V.emplace_back(x + s, y + s);
  V.emplace_back(x, y + s);
  if (s > 1) {
    const auto mask = qt.edge_mask_smaller(L);
    if (mask.test(0))
      V.emplace_back(x + (s >> 1), y);
    if (mask.test(1))
      V.emplace_back(x + s, y + (s >> 1));
    if (mask.test(2))
      V.emplace_back(x + (s >> 1), y + s);
    if (mask.test(3))
      V.emplace_back(x, y + (s >> 1));
  }
  return V;
}

// Check that all triangle vertices belong to the allowed boundary set of the
// given leaf
template <class TriangulationT>
static bool triangle_vertices_within_leaf_boundary(
  const NonUniformQuadtreeLocator &qt,
  const NonUniformQuadtreeLocator::Leaf &L,
  const typename TriangulationT::Triangle &T)
{
  const auto V = leaf_allowed_vertices(qt, L);
  auto is_allowed = [&](uint32_t vx, uint32_t vy) {
    for (auto [ax, ay] : V)
      if (vx == ax && vy == ay)
        return true;
    return false;
  };
  for (size_t i = 0; i < 3; ++i) {
    const uint32_t vx = T.vertices[i].x();
    const uint32_t vy = T.vertices[i].y();
    if (!is_allowed(vx, vy))
      return false;
  }
  return true;
}

static std::vector<Point> sample_uniform_points(size_t N)
{
  std::vector<Point> out;
  out.reserve(N);
  std::mt19937_64 rng(0xD00DFEEDULL);
  std::uniform_real_distribution<double> U(1e-7, 8.0 - 1e-7);
  for (size_t i = 0; i < N; ++i)
    out.emplace_back(Point(U(rng), U(rng)));
  return out;
}

// Adversarial samples near hanging nodes (on both sides), plus leaf centers
static std::vector<Point> sample_critical_points(
  const NonUniformQuadtreeLocator &qt)
{
  std::vector<Point> pts;
  const auto &Ls = qt.leaves();
  for (const auto &L : Ls) {
    const double x = L.x, y = L.y, s = L.size;

    // center
    pts.emplace_back(Point(x + 0.5 * s, y + 0.5 * s));

    if (s <= 1)
      continue;
    const double eps = 1e-6;

    const auto mask = qt.edge_mask_smaller(L);
    if (mask.test(0)) {  // bottom
      if (y > 0) {
        pts.emplace_back(Point(x + 0.5 * s, y + eps));  // just inside big leaf
        pts.emplace_back(
          Point(x + 0.25 * s, y - eps));  // just inside the smaller neighbor
      }
    }
    if (mask.test(1)) {  // right
      if (x + s < GRID_W) {
        pts.emplace_back(Point(x + s - eps, y + 0.5 * s));
        pts.emplace_back(Point(x + s + eps, y + 0.25 * s));
      }
    }
    if (mask.test(2)) {  // top
      if (y + s < GRID_H) {
        pts.emplace_back(Point(x + 0.5 * s, y + s - eps));
        pts.emplace_back(Point(x + 0.25 * s, y + s + eps));
      }
    }
    if (mask.test(3)) {  // left
      if (x > 0) {
        pts.emplace_back(Point(x + eps, y + 0.5 * s));
        pts.emplace_back(Point(x - eps, y + 0.25 * s));
      }
    }
  }
  return pts;
}

template <class Proj>
static void run_case(
  const char *label,
  const NonUniformQuadtreeLocator &qt,
  const Proj &proj,
  bool expect_ok)
{
  BOOST_TEST_CONTEXT(label)
  {
    Triangulation<NonUniformQuadtreeLocator, Proj> tri;
    const bool ok = tri.build(&qt, &proj);
    BOOST_CHECK_EQUAL(ok, expect_ok);
    if (!ok)
      return;

    const auto &Ts = tri.triangles();

    // Expected triangle count: sum over leaves of (n-2), with n = 4 +
    // (#hanging nodes)
    size_t expected = 0;
    for (const auto &L : qt.leaves())
      expected += static_cast<size_t>(2 + qt.count_midpoints(L));
    BOOST_CHECK_EQUAL(Ts.size(), expected);

    // Every triangle's vertices must belong to its owning leaf boundary set
    // (we identify the owner by the triangle centroid -> quadtree locate)
    using TriT = Triangulation<NonUniformQuadtreeLocator, Proj>;
    for (const auto &T : Ts) {
      const double cx =
        (T.vertices[0].x() + T.vertices[1].x() + T.vertices[2].x()) / 3.0;
      const double cy =
        (T.vertices[0].y() + T.vertices[1].y() + T.vertices[2].y()) / 3.0;
      const auto L = qt.locate(Point(cx, cy));
      BOOST_CHECK(triangle_vertices_within_leaf_boundary<TriT>(qt, L, T));
    }

    // Locate() correctness on deterministic sets
    auto try_points = [&](const std::vector<Point> &P) {
      for (const auto &p : P) {
        auto it = tri.locate(p);
        BOOST_REQUIRE(it != Ts.end());

        // Containment under EPICK logic
        BOOST_CHECK(triangle_contains<TriT>(*it, p));

        // Triangle should belong to the same leaf as quadtree locate(p)
        const auto leaf = qt.locate(p);
        BOOST_CHECK(
          triangle_vertices_within_leaf_boundary<TriT>(qt, leaf, *it));

        // Cache fast-path: re-locate same point -> same index
        const auto idx1 = static_cast<uint32_t>(std::distance(Ts.begin(), it));
        const auto it2 = tri.locate(p);
        const auto idx2 =
          static_cast<uint32_t>(std::distance(Ts.begin(), it2));
        BOOST_CHECK_EQUAL(idx1, idx2);
      }
    };

    try_points(sample_uniform_points(2000));
    try_points(sample_critical_points(qt));
  }
}

BOOST_AUTO_TEST_CASE(NonUniform_A_Shear_Positive)
{
  auto qt = NonUniformQuadtreeLocator::make_layout_A();
  BalancedProjection<ShearX> proj(&qt, ShearX{+0.15});
  run_case("A: Shear +X k=0.15", qt, proj, /*expect_ok=*/true);
}

BOOST_AUTO_TEST_CASE(NonUniform_A_Shear_Negative)
{
  auto qt = NonUniformQuadtreeLocator::make_layout_A();
  BalancedProjection<ShearX> proj(&qt, ShearX{-0.22});
  run_case("A: Shear +X k=-0.22", qt, proj, /*expect_ok=*/true);
}

BOOST_AUTO_TEST_CASE(NonUniform_A_Rotate)
{
  auto qt = NonUniformQuadtreeLocator::make_layout_A();
  BalancedProjection<Rotate> proj(&qt, Rotate{M_PI / 9.0});  // ~20 deg
  run_case("A: Rotate 20deg", qt, proj, /*expect_ok=*/true);
}

BOOST_AUTO_TEST_CASE(NonUniform_A_Radial_Shrink)
{
  auto qt = NonUniformQuadtreeLocator::make_layout_A();
  BalancedProjection<RadialInward> proj(&qt, RadialInward{0.30});
  run_case("A: Radial inward alpha=0.30", qt, proj, /*expect_ok=*/true);
}

BOOST_AUTO_TEST_CASE(NonUniform_A_Wavy)
{
  auto qt = NonUniformQuadtreeLocator::make_layout_A();
  BalancedProjection<Wavy> proj(&qt, Wavy{0.05});
  run_case("A: Wavy eps=0.05", qt, proj, /*expect_ok=*/true);
}

BOOST_AUTO_TEST_CASE(NonUniform_A_Flip_Is_Detected)
{
  auto qt = NonUniformQuadtreeLocator::make_layout_A();
  BalancedProjection<MirrorXAdapter> proj(&qt, MirrorXAdapter{});
  run_case("A: Mirror X (flip expected)", qt, proj, /*expect_ok=*/false);
}

BOOST_AUTO_TEST_CASE(NonUniform_B_Wavy)
{
  auto qt = NonUniformQuadtreeLocator::make_layout_B();
  BalancedProjection<Wavy> proj(&qt, Wavy{0.07});
  run_case("B: Wavy eps=0.07", qt, proj, /*expect_ok=*/true);
}

BOOST_AUTO_TEST_CASE(NonUniform_B_Rotate_And_Shear)
{
  auto qt = NonUniformQuadtreeLocator::make_layout_B();
  // compose mild rotate + shear by nesting in a lambda-based deform (kept
  // small and orientation-preserving)
  struct Combo {
    std::pair<double, double> operator()(
      double x,
      double y,
      double w,
      double h) const
    {
      // rotate
      const double cx = 0.5 * w, cy = 0.5 * h, c = std::cos(M_PI / 12.0),
                   s = std::sin(M_PI / 12.0);
      const double X = x - cx, Y = y - cy;
      double xr = c * X - s * Y + cx, yr = s * X + c * Y + cy;
      // shear
      xr += 0.1 * yr;
      return {xr, yr};
    }
  };
  BalancedProjection<Combo> proj(&qt, Combo{});
  run_case("B: Rotate 15deg + Shear 0.1", qt, proj, /*expect_ok=*/true);
}
