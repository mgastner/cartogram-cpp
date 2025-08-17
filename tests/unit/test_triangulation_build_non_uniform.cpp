#define BOOST_TEST_MODULE test_triangulation_build_non_uniform
#include "cgal_typedef.hpp"
#include "quadtree_corner.hpp"
#include "round_point.hpp"
#include "triangulation.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <algorithm>
#include <array>
#include <bitset>
#include <boost/test/included/unit_test.hpp>
#include <cstdint>
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
  };  // bottom-left, size in {1,2,4}

  NonUniformQuadtreeLocator() = default;

  static NonUniformQuadtreeLocator make_layout_A()
  {
    // Mix of 4x4, 2x2, and a 1x1 patch
    std::vector<Leaf> L;

    // TR 4x4
    L.push_back({4, 4, 4});

    // TL four 2x2
    L.push_back({0, 4, 2});
    L.push_back({2, 4, 2});
    L.push_back({0, 6, 2});
    L.push_back({2, 6, 2});

    // BR four 2x2
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

  // Bitset edges with smaller neighbor (0: bottom, 1: right, 2: top, 3: left)
  std::bitset<4> edge_mask_smaller(const Leaf &L) const
  {
    std::bitset<4> m;
    const uint32_t x = L.x, y = L.y, s = L.size;
    if (s > 1) {
      if (y > 0)
        m.set(0, size_at_(x, y - 1) < s);
      if (x + s < w_)
        m.set(1, size_at_(x + s, y) < s);
      if (y + s < h_)
        m.set(2, size_at_(x, y + s) < s);
      if (x > 0)
        m.set(3, size_at_(x - 1, y) < s);
    }
    return m;
  }

private:
  uint32_t w_{GRID_W}, h_{GRID_H};
  std::vector<Leaf> leaves_;
  int cell_leaf_idx_[GRID_W][GRID_H]{};
  uint32_t cell_leaf_size_[GRID_W][GRID_H]{};

  void build_maps_()
  {
    for (uint32_t y = 0; y < h_; ++y)
      for (uint32_t x = 0; x < w_; ++x)
        cell_leaf_idx_[x][y] = -1, cell_leaf_size_[x][y] = 0;

    for (size_t i = 0; i < leaves_.size(); ++i) {
      const auto &L = leaves_[i];
      for (uint32_t yy = L.y; yy < L.y + L.size; ++yy)
        for (uint32_t xx = L.x; xx < L.x + L.size; ++xx) {
          BOOST_REQUIRE(cell_leaf_idx_[xx][yy] == -1);
          cell_leaf_idx_[xx][yy] = static_cast<int>(i);
          cell_leaf_size_[xx][yy] = L.size;
        }
    }
    for (uint32_t y = 0; y < h_; ++y)
      for (uint32_t x = 0; x < w_; ++x)
        BOOST_REQUIRE(cell_leaf_idx_[x][y] >= 0);
  }

  uint32_t size_at_(uint32_t cx, uint32_t cy) const
  {
    return cell_leaf_size_[cx][cy];
  }
};

// pack (x,y) for hash keys
static inline uint64_t pack_xy(uint32_t x, uint32_t y)
{
  return (uint64_t(x) << 32) | uint64_t(y);
}

// key64 for (leaf_x, leaf_y, leaf_size) -> 64-bit
static inline uint64_t key64(uint32_t x, uint32_t y, uint32_t s)
{
  return (uint64_t(x) << 40) | (uint64_t(y) << 20) | uint64_t(s);
}

// Deform: Shear in X: x' = x + k*y
struct ShearX {
  double k{};
};

template <class Deform> struct BalancedProjection {
  BalancedProjection(const NonUniformQuadtreeLocator *qt, Deform deform)
      : qt_(qt), deform_(std::move(deform)), w_(qt->width()), h_(qt->height())
  {
    const auto &leaves = qt_->leaves();

    for (size_t i = 0; i < leaves.size(); ++i) {
      const auto &L = leaves[i];
      offset_map_[pack_xy(L.x, L.y)] = static_cast<uint32_t>(i * 4);
    }

    // Mark midpoints on edges where neighbor is smaller (+1 higher depth)
    for (const auto &L : leaves) {
      const uint32_t x = L.x, y = L.y, s = L.size;
      if (s <= 1)
        continue;
      const auto mask = qt_->edge_mask_smaller(L);
      if (mask.test(0))
        mid_.insert(pack_xy(x + (s >> 1), y));
      if (mask.test(1))
        mid_.insert(pack_xy(x + s, y + (s >> 1)));
      if (mask.test(2))
        mid_.insert(pack_xy(x + (s >> 1), y + s));
      if (mask.test(3))
        mid_.insert(pack_xy(x, y + (s >> 1)));
    }
  }

  bool is_valid_corner(uint32_t x, uint32_t y) const
  {
    return mid_.find(pack_xy(x, y)) != mid_.end();
  }

  uint32_t offset(uint32_t x, uint32_t y) const
  {
    auto it = offset_map_.find(pack_xy(x, y));
    BOOST_REQUIRE(it != offset_map_.end());
    return it->second;
  }

  Point get(uint32_t x, uint32_t y) const
  {
    const double xd = static_cast<double>(x);
    const double yd = static_cast<double>(y);
    const double k = deform_.k;
    return Point(xd + k * yd, yd);
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
  std::unordered_set<uint64_t> mid_;
};

// Layout A (shear +0.15): mixture of 4x4, 2x2, 1x1 (balanced).
// - TR: one 4x4 at (4,4)
// - TL: four 2x2 over [0,4)x[4,8)
// - BR: four 2x2 over [4,8)x[0,4)
// - BL: mostly 2x2 with a 2x2 block of 1x1 at (0..2,0..2)
static const std::unordered_map<uint64_t, std::vector<std::array<int, 3>>>
  GOLD_A_SHEAR_POS = {
    {0x0000000000000001ULL, {{0, 1, 3}, {1, 2, 3}}},
    {0x0000000000100001ULL, {{0, 1, 3}, {1, 2, 3}}},
    {0x0000000000200002ULL, {{0, 1, 4}, {1, 3, 4}, {1, 2, 3}}},
    {0x0000000000400002ULL, {{0, 1, 3}, {1, 2, 3}}},
    {0x0000000000600002ULL, {{0, 1, 3}, {1, 2, 3}}},
    {0x0000010000000001ULL, {{0, 1, 3}, {1, 2, 3}}},
    {0x0000010000100001ULL, {{0, 1, 3}, {1, 2, 3}}},
    {0x0000020000000002ULL, {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}}},
    {0x0000020000200002ULL, {{0, 1, 3}, {1, 2, 3}}},
    {0x0000020000400002ULL, {{0, 1, 3}, {1, 2, 3}}},
    {0x0000020000600002ULL, {{0, 1, 3}, {1, 2, 3}}},
    {0x0000040000000002ULL, {{0, 1, 3}, {1, 2, 3}}},
    {0x0000040000200002ULL, {{0, 1, 3}, {1, 2, 3}}},
    {0x0000040000400004ULL, {{0, 1, 5}, {1, 3, 5}, {1, 2, 3}, {3, 4, 5}}},
    {0x0000060000000002ULL, {{0, 1, 3}, {1, 2, 3}}},
    {0x0000060000200002ULL, {{0, 1, 3}, {1, 2, 3}}},
};

template <class Proj>
static void build_leaf_polygon(
  const NonUniformQuadtreeLocator::Leaf &L,
  const Proj &proj,
  std::vector<QuadtreeCorner> &poly)
{
  const uint32_t x = L.x, y = L.y, s = L.size;
  poly.clear();
  poly.reserve(8);
  auto push = [&](uint32_t xx, uint32_t yy) {
    poly.push_back({xx, yy});
  };

  push(x, y);
  if (s > 1 && proj.is_valid_corner(x + (s >> 1), y))
    push(x + (s >> 1), y);
  push(x + s, y);
  if (s > 1 && proj.is_valid_corner(x + s, y + (s >> 1)))
    push(x + s, y + (s >> 1));
  push(x + s, y + s);
  if (s > 1 && proj.is_valid_corner(x + (s >> 1), y + s))
    push(x + (s >> 1), y + s);
  push(x, y + s);
  if (s > 1 && proj.is_valid_corner(x, y + (s >> 1)))
    push(x, y + (s >> 1));
}

// Map a triangle's 3 integer vertices to indices in `poly`
static std::array<int, 3> tri_to_poly_indices(
  const std::vector<QuadtreeCorner> &poly,
  const typename Triangulation<
    NonUniformQuadtreeLocator,
    BalancedProjection<ShearX>>::Triangle &T)
{
  auto idx_of = [&](uint32_t vx, uint32_t vy) -> int {
    for (size_t i = 0; i < poly.size(); ++i)
      if (poly[i].x() == vx && poly[i].y() == vy)
        return static_cast<int>(i);
    return -1;
  };

  std::array<int, 3> out{};
  for (size_t i = 0; i < 3; ++i) {
    out[i] = idx_of(T.vertices[i].x(), T.vertices[i].y());
  }
  return out;
}

static std::array<int, 3> canon(std::array<int, 3> t)
{
  std::sort(t.begin(), t.end());
  return t;
}

template <class Proj>
static void run_gold_case(
  const char *label,
  const NonUniformQuadtreeLocator &qt,
  const Proj &proj,
  const std::unordered_map<uint64_t, std::vector<std::array<int, 3>>> &gold)
{
  BOOST_TEST_CONTEXT(label)
  {
    Triangulation<NonUniformQuadtreeLocator, Proj> tri;
    const bool ok = tri.build(&qt, &proj);
    BOOST_REQUIRE(ok);

    const auto &Ts = tri.triangles();

    // Expected global triangle count = sum of golden per-leaf counts
    size_t expected_total = 0;
    for (auto &kv : gold)
      expected_total += kv.second.size();
    BOOST_CHECK_EQUAL(Ts.size(), expected_total);

    // For each leaf, compare its triangle set against the golden set
    std::vector<QuadtreeCorner> poly;
    for (const auto &L : qt.leaves()) {
      const uint64_t key = key64(L.x, L.y, L.size);
      auto itg = gold.find(key);
      BOOST_REQUIRE_MESSAGE(
        itg != gold.end(),
        "Missing golden entry for leaf (" << L.x << "," << L.y
                                          << "; s=" << L.size << ")");

      build_leaf_polygon(L, proj, poly);

      // Gather triangles that belong to this leaf (via centroid locate)
      std::vector<std::array<int, 3>> got_idxs;
      for (const auto &T : Ts) {
        const double cx =
          (T.vertices[0].x() + T.vertices[1].x() + T.vertices[2].x()) / 3.0;
        const double cy =
          (T.vertices[0].y() + T.vertices[1].y() + T.vertices[2].y()) / 3.0;
        const auto Lhat = qt.locate(Point(cx, cy));
        if (!(Lhat.x == L.x && Lhat.y == L.y && Lhat.size == L.size))
          continue;

        auto tri_idx = tri_to_poly_indices(poly, T);
        BOOST_REQUIRE(tri_idx[0] >= 0 && tri_idx[1] >= 0 && tri_idx[2] >= 0);
        got_idxs.push_back(canon(tri_idx));
      }

      // Canonicalize expected as well
      std::vector<std::array<int, 3>> exp_idxs;
      exp_idxs.reserve(itg->second.size());
      for (auto t : itg->second)
        exp_idxs.push_back(canon(t));

      auto sort3 =
        [](const std::array<int, 3> &a, const std::array<int, 3> &b) {
          if (a[0] != b[0])
            return a[0] < b[0];
          if (a[1] != b[1])
            return a[1] < b[1];
          return a[2] < b[2];
        };
      std::sort(got_idxs.begin(), got_idxs.end(), sort3);
      std::sort(exp_idxs.begin(), exp_idxs.end(), sort3);

      // clang-format off
      BOOST_CHECK_EQUAL(got_idxs.size(), exp_idxs.size());
      BOOST_CHECK(
        std::is_permutation(
          got_idxs.begin(),
          got_idxs.end(),
          exp_idxs.begin(),
          exp_idxs.end()));
      // clang-format on
    }
  }
}

BOOST_AUTO_TEST_CASE(LayoutA_Shear_Positive_HardcodedOptimal)
{
  auto qt = NonUniformQuadtreeLocator::make_layout_A();
  BalancedProjection<ShearX> proj(&qt, ShearX{+0.15});
  run_gold_case("A: Shear +X k=0.15 hardcoded", qt, proj, GOLD_A_SHEAR_POS);
}

// TODO: Add more test cases for different layouts and projections
