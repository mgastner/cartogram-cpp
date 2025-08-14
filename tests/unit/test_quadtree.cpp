#define BOOST_TEST_MODULE test_quadtree
#include "quadtree.hpp"
#include "quadtree_leaf_locator.hpp"
#include <bit>
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <cstdlib>
#include <vector>

BOOST_AUTO_TEST_SUITE(QuadtreeTests)

BOOST_AUTO_TEST_CASE(CTAD_and_initial_state)
{
  auto metric = [](uint32_t, uint32_t, uint32_t s) {
    return s;
  };
  Quadtree qt(16, 1, metric);
  BOOST_TEST(qt.num_leaves() == 1u);
}

BOOST_AUTO_TEST_CASE(Build_reaches_target_within_bounds)
{
  constexpr std::size_t target = 50;
  auto metric = [](uint32_t, uint32_t, uint32_t s) {
    return s;
  };

  Quadtree qt(256, target, metric);
  qt.build();

  const std::size_t n = qt.num_leaves();
  BOOST_TEST(n >= target);
  BOOST_TEST(n <= target + 3);
}

BOOST_AUTO_TEST_CASE(Build_does_not_oversplit_unit_cells)
{
  auto metric = [](uint32_t, uint32_t, uint32_t) {
    return 1.0f;
  };
  Quadtree qt(4, 100, metric);  // 4 * 4 domain -> 16 unit cells max
  qt.build();
  BOOST_TEST(qt.num_leaves() == 16u);
}

BOOST_AUTO_TEST_CASE(Leaves_form_exact_tiling)
{
  constexpr uint32_t root = 16;
  auto metric = [](uint32_t, uint32_t, uint32_t s) {
    return s;
  };

  Quadtree qt(root, 40, metric);
  qt.build();

  std::vector<uint8_t> mask(root * root, 0);
  for (auto [x, y, sz] : qt.leaves()) {
    for (uint32_t yy = y; yy < y + sz; ++yy)
      for (uint32_t xx = x; xx < x + sz; ++xx) {
        BOOST_TEST(mask[yy * root + xx] == 0u);  // no overlap
        mask[yy * root + xx] = 1;
      }
  }
  for (uint8_t v : mask)
    BOOST_TEST(v == 1u);  // full coverage
}

BOOST_AUTO_TEST_CASE(Grading_keeps_depth_difference_within_one)
{
  constexpr uint32_t root = 128;
  auto metric = [](uint32_t x, uint32_t y, uint32_t s) {
    return double((x ^ y) + s);  // mildly irregular priority
  };

  Quadtree qt(root, 300, metric);
  qt.build();
  qt.grade();

  const auto leaves = qt.leaves();
  auto depth_of = [=](uint32_t sz) {
    return int(std::countr_zero(root) - std::countr_zero(sz));
  };

  const std::size_t n = leaves.size();
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      const auto a = leaves[i];
      const auto b = leaves[j];

      const bool overlapY = a.y < b.y + b.size && b.y < a.y + a.size;
      const bool overlapX = a.x < b.x + b.size && b.x < a.x + a.size;

      const bool adjacent =
        ((a.x + a.size == b.x || b.x + b.size == a.x) && overlapY) ||
        ((a.y + a.size == b.y || b.y + b.size == a.y) && overlapX);

      if (adjacent) {
        const int da = depth_of(a.size);
        const int db = depth_of(b.size);
        BOOST_TEST(std::abs(da - db) <= 1);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(Locate_clamps_and_boundaries)
{
  constexpr uint32_t root = 16;
  auto metric = [](uint32_t, uint32_t, uint32_t) {
    return 1.0;
  };
  Quadtree qt(root, 10000, metric);  // force full refinement
  qt.build();

  QuadtreeLeafLocator qt_locator;
  qt_locator.build(qt.root_size(), qt.nodes());

  auto a = qt_locator.locate(-10.0, -5.0);
  BOOST_TEST(a.x == 0u);
  BOOST_TEST(a.y == 0u);
  BOOST_TEST(a.size == 1u);

  auto b = qt_locator.locate(root, 0.0);
  BOOST_TEST(b.x == root - 1);
  BOOST_TEST(b.y == 0u);
  BOOST_TEST(b.size == 1u);

  auto c = qt_locator.locate(0.0, root);
  BOOST_TEST(c.x == 0u);
  BOOST_TEST(c.y == root - 1);
  BOOST_TEST(c.size == 1u);

  auto d = qt_locator.locate(root, root);
  BOOST_TEST(d.x == root - 1);
  BOOST_TEST(d.y == root - 1);
  BOOST_TEST(d.size == 1u);

  auto e = qt_locator.locate(8.0, 8.0);
  BOOST_TEST(e.x == 8u);
  BOOST_TEST(e.y == 8u);
  BOOST_TEST(e.size == 1u);
}

BOOST_AUTO_TEST_CASE(Graded_leaves_form_exact_tiling_and_alignment)
{
  constexpr uint32_t root = 128;
  auto metric = [](uint32_t, uint32_t, uint32_t s) {
    return s;
  };

  Quadtree qt(root, 300, metric);
  qt.build();
  qt.grade();

  std::vector<uint8_t> mask(root * root, 0);
  for (auto [x, y, sz] : qt.leaves()) {
    BOOST_TEST(std::has_single_bit(sz));
    BOOST_TEST(x % sz == 0u);
    BOOST_TEST(y % sz == 0u);

    for (uint32_t yy = y; yy < y + sz; ++yy)
      for (uint32_t xx = x; xx < x + sz; ++xx) {
        BOOST_TEST(mask[yy * root + xx] == 0u);  // no overlap
        mask[yy * root + xx] = 1;
      }
  }
  for (uint8_t v : mask)
    BOOST_TEST(v == 1u);  // full coverage
}

BOOST_AUTO_TEST_CASE(Locate_covers_all_pixels_after_build_and_grade)
{
  constexpr uint32_t root = 64;
  auto metric = [](uint32_t x, uint32_t y, uint32_t s) {
    return double((x ^ y) + s);  // irregular
  };
  Quadtree qt(root, 500, metric);
  qt.build();
  qt.grade();

  QuadtreeLeafLocator qt_locator;
  qt_locator.build(qt.root_size(), qt.nodes());

  for (uint32_t py = 0; py < root; ++py) {
    for (uint32_t px = 0; px < root; ++px) {
      auto lf = qt_locator.locate(double(px), double(py));
      BOOST_TEST(px >= lf.x);
      BOOST_TEST(px < lf.x + lf.size);
      BOOST_TEST(py >= lf.y);
      BOOST_TEST(py < lf.y + lf.size);
    }
  }
}

BOOST_AUTO_TEST_CASE(Priority_driven_refinement)
{
  // Integral metric; strong pull to (0,0) only
  constexpr uint32_t root = 64;
  auto metric = [](uint32_t x, uint32_t y, uint32_t) -> uint32_t {
    return (x == 0u && y == 0u) ? 1u : 0u;
  };

  Quadtree qt(root, 100, metric);
  qt.build();

  QuadtreeLeafLocator qt_locator;
  qt_locator.build(qt.root_size(), qt.nodes());

  auto a = qt_locator.locate(0.0, 0.0);
  BOOST_TEST(a.x == 0u);
  BOOST_TEST(a.y == 0u);
  BOOST_TEST(a.size == 1u);  // refined deeply at the hotspot

  auto b = qt_locator.locate(root - 1.0, root - 1.0);
  BOOST_TEST(b.size > 1u);  // away from hotspot stays larger
}

BOOST_AUTO_TEST_CASE(Grade_balances_neighbors_next_to_unit_cell)
{
  constexpr uint32_t root = 64;
  auto metric = [](uint32_t x, uint32_t y, uint32_t) -> uint32_t {
    return (x == 0u && y == 0u) ? 1u : 0u;
  };

  Quadtree qt(root, 100, metric);
  qt.build();

  QuadtreeLeafLocator qt_locator;
  qt_locator.build(qt.root_size(), qt.nodes());

  // Ensure (0,0) is unit-sized before grading
  auto a = qt_locator.locate(0.0, 0.0);
  BOOST_TEST(a.size == 1u);

  qt.grade();

  // Neighbors along +x and +y must become size <= 2 for 1‑irregularity
  auto nx = qt_locator.locate(1.0, 0.0);
  auto ny = qt_locator.locate(0.0, 1.0);
  BOOST_TEST(nx.size <= 2u);
  BOOST_TEST(ny.size <= 2u);
}

BOOST_AUTO_TEST_CASE(NumLeaves_consistent_with_leaves_vector)
{
  constexpr uint32_t root = 32;
  auto metric = [](uint32_t, uint32_t, uint32_t s) {
    return s;
  };
  Quadtree qt(root, 200, metric);
  qt.build();
  BOOST_TEST(qt.leaves().size() == qt.num_leaves());
  qt.grade();
  BOOST_TEST(qt.leaves().size() == qt.num_leaves());
}

BOOST_AUTO_TEST_CASE(Locate_consistent_within_leaf)
{
  constexpr uint32_t root = 64;
  auto metric = [](uint32_t, uint32_t, uint32_t s) {
    return s;
  };
  Quadtree qt(root, 150, metric);
  qt.build();

  QuadtreeLeafLocator qt_locator;
  qt_locator.build(qt.root_size(), qt.nodes());

  auto lf = qt_locator.locate(20.2, 37.7);
  const double eps = 1e-9;
  auto p1 = qt_locator.locate(lf.x + 0.1, lf.y + 0.1);
  auto p2 = qt_locator.locate(lf.x + lf.size - eps, lf.y + 0.1);
  auto p3 = qt_locator.locate(lf.x + 0.1, lf.y + lf.size - eps);
  auto p4 = qt_locator.locate(lf.x + lf.size - eps, lf.y + lf.size - eps);

  auto same = [](auto a, auto b) {
    return a.x == b.x && a.y == b.y && a.size == b.size;
  };
  BOOST_TEST(same(p1, lf));
  BOOST_TEST(same(p2, lf));
  BOOST_TEST(same(p3, lf));
  BOOST_TEST(same(p4, lf));
}

BOOST_AUTO_TEST_CASE(Grade_noop_on_fully_refined_uniform_tree)
{
  constexpr uint32_t root = 8;
  auto metric = [](uint32_t, uint32_t, uint32_t) {
    return 1.0;
  };
  Quadtree qt(root, 1000, metric);  // enough to reach unit cells
  qt.build();

  const auto before = qt.leaves();
  const std::size_t n = qt.num_leaves();

  qt.grade();

  const auto after = qt.leaves();
  BOOST_TEST(qt.num_leaves() == n);
  BOOST_TEST(before.size() == after.size());
  for (std::size_t i = 0; i < before.size(); ++i) {
    BOOST_TEST(before[i].x == after[i].x);
    BOOST_TEST(before[i].y == after[i].y);
    BOOST_TEST(before[i].size == after[i].size);
  }
}

BOOST_AUTO_TEST_CASE(Predictable_4x4_build_and_grade)
{
  constexpr uint32_t root = 4;

  // After the root split, only the top-left 2x2 gets split
  // again
  auto metric = [](uint32_t x, uint32_t y, uint32_t s) -> int {
    if (s == 2 && x == 0u && y == 0u)
      return 2;  // prefer (0,0) child
    if (s == 2)
      return 1;  // other 2x2 children
    return 0;  // everything else
  };

  using QT = Quadtree<decltype(metric)>;

  QT qt(
    root,
    /*target*/ 7,
    metric);  // 1 -> split root to 4 -> split one child -> 7 leaves total
  qt.build();

  // Expected leaves after build:
  //   - TL quadrant fully refined: (0,0),(1,0),(0,1),(1,1) size=1
  //   - Other quadrants: size=2 at (2,0), (0,2), (2,2)
  std::vector<QT::Leaf> expected{
    {0, 0, 1},
    {1, 0, 1},
    {0, 1, 1},
    {1, 1, 1},
    {2, 0, 2},
    {0, 2, 2},
    {2, 2, 2}};

  auto actual = qt.leaves();
  auto cmp = [](const QT::Leaf &a, const QT::Leaf &b) {
    if (a.y != b.y)
      return a.y < b.y;
    if (a.x != b.x)
      return a.x < b.x;
    return a.size < b.size;
  };
  std::sort(expected.begin(), expected.end(), cmp);
  std::sort(actual.begin(), actual.end(), cmp);

  BOOST_TEST(actual.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    BOOST_TEST(actual[i].x == expected[i].x);
    BOOST_TEST(actual[i].y == expected[i].y);
    BOOST_TEST(actual[i].size == expected[i].size);
  }

  // Manual tiling of all 16 pixels
  std::vector<unsigned char> mask(root * root, 0);
  for (const auto &lf : actual) {
    for (uint32_t yy = lf.y; yy < lf.y + lf.size; ++yy)
      for (uint32_t xx = lf.x; xx < lf.x + lf.size; ++xx)
        ++mask[yy * root + xx];
  }
  for (unsigned char v : mask)
    BOOST_TEST(v == 1);

  // Check locate() at pixel centers for all 16 pixels
  auto expected_for = [](uint32_t px, uint32_t py) -> QT::Leaf {
    if (px < 2 && py < 2)
      return QT::Leaf{px, py, 1};
    if (px >= 2 && py < 2)
      return QT::Leaf{2, 0, 2};
    if (px < 2 && py >= 2)
      return QT::Leaf{0, 2, 2};
    return QT::Leaf{2, 2, 2};
  };

  QuadtreeLeafLocator qt_locator;
  qt_locator.build(qt.root_size(), qt.nodes());

  for (uint32_t py = 0; py < root; ++py) {
    for (uint32_t px = 0; px < root; ++px) {
      auto lf = qt_locator.locate(px + 0.5, py + 0.5);
      auto ex = expected_for(px, py);
      BOOST_TEST(lf.x == ex.x);
      BOOST_TEST(lf.y == ex.y);
      BOOST_TEST(lf.size == ex.size);
    }
  }

  // grade() is a no-op (configuration already 1‑irregular); re-check
  qt.grade();

  auto after = qt.leaves();
  std::sort(after.begin(), after.end(), cmp);
  BOOST_TEST(after.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    BOOST_TEST(after[i].x == expected[i].x);
    BOOST_TEST(after[i].y == expected[i].y);
    BOOST_TEST(after[i].size == expected[i].size);
  }

  qt_locator.build(qt.root_size(), qt.nodes());

  for (uint32_t py = 0; py < root; ++py) {
    for (uint32_t px = 0; px < root; ++px) {
      auto lf = qt_locator.locate(px + 0.5, py + 0.5);
      auto ex = expected_for(px, py);
      BOOST_TEST(lf.x == ex.x);
      BOOST_TEST(lf.y == ex.y);
      BOOST_TEST(lf.size == ex.size);
    }
  }
}

BOOST_AUTO_TEST_CASE(Predictable_8x8_build_then_grade_changes)
{
  constexpr uint32_t root = 8;

  // Build plan:
  // 1) Split root (4x4 children). Prefer TL (0,0,4)
  // 2) Split TL's 2x2 at (2,0) so we get four 1x1 touching the TR 4x4
  // boundary
  auto metric = [](uint32_t x, uint32_t y, uint32_t s) -> int {
    if (s == 4)
      return (x == 0u && y == 0u) ? 100 : 1;  // choose TL 4x4 first
    if (s == 2)
      return (x == 2u && y == 0u) ? 1000 : 1;  // then pick (2,0) 2x2
    return 0;
  };

  Quadtree qt(root, /*target*/ 10, metric);  // 1->4->7->10 leaves
  qt.build();

  auto leaves = qt.leaves();
  BOOST_TEST(leaves.size() == 10u);

  auto has = [&](uint32_t x, uint32_t y, uint32_t sz) {
    for (auto &lf : leaves)
      if (lf.x == x && lf.y == y && lf.size == sz)
        return true;
    return false;
  };

  // After build:
  // - 4x4 at TR, BL, BR
  BOOST_TEST(has(4, 0, 4));
  BOOST_TEST(has(0, 4, 4));
  BOOST_TEST(has(4, 4, 4));
  // - TL has 2x2 except (2,0) which was split
  BOOST_TEST(has(0, 0, 2));
  BOOST_TEST(has(0, 2, 2));
  BOOST_TEST(has(2, 2, 2));
  // - 1x1 from splitting (2,0)
  BOOST_TEST(has(2, 0, 1));
  BOOST_TEST(has(3, 0, 1));
  BOOST_TEST(has(2, 1, 1));
  BOOST_TEST(has(3, 1, 1));

  // Full coverage + locate() check
  std::vector<unsigned char> mask(root * root, 0);
  for (auto lf : leaves)
    for (uint32_t yy = lf.y; yy < lf.y + lf.size; ++yy)
      for (uint32_t xx = lf.x; xx < lf.x + lf.size; ++xx)
        ++mask[yy * root + xx];
  for (unsigned char v : mask)
    BOOST_TEST(v == 1);

  QuadtreeLeafLocator qt_locator;
  qt_locator.build(qt.root_size(), qt.nodes());
  for (uint32_t py = 0; py < root; ++py)
    for (uint32_t px = 0; px < root; ++px) {
      auto lf = qt_locator.locate(px + 0.5, py + 0.5);
      BOOST_TEST(px >= lf.x);
      BOOST_TEST(px < lf.x + lf.size);
      BOOST_TEST(py >= lf.y);
      BOOST_TEST(py < lf.y + lf.size);
    }

  // Grade: 1x1 next to TR 4x4 -> depth diff = 2 -> TR must split into four 2x2
  qt.grade();

  auto after = qt.leaves();
  BOOST_TEST(after.size() == 13u);

  auto has_after = [&](uint32_t x, uint32_t y, uint32_t sz) {
    for (auto &lf : after)
      if (lf.x == x && lf.y == y && lf.size == sz)
        return true;
    return false;
  };

  // TR 4x4 gone, replaced by four 2x2
  BOOST_TEST(!has_after(4, 0, 4));
  BOOST_TEST(has_after(4, 0, 2));
  BOOST_TEST(has_after(6, 0, 2));
  BOOST_TEST(has_after(4, 2, 2));
  BOOST_TEST(has_after(6, 2, 2));

  // Spot‑check survivors
  BOOST_TEST(has_after(0, 0, 2));
  BOOST_TEST(has_after(3, 1, 1));
  BOOST_TEST(has_after(4, 4, 4));

  // Coverage + locate again
  std::fill(mask.begin(), mask.end(), 0);
  for (auto lf : after)
    for (uint32_t yy = lf.y; yy < lf.y + lf.size; ++yy)
      for (uint32_t xx = lf.x; xx < lf.x + lf.size; ++xx)
        ++mask[yy * root + xx];
  for (unsigned char v : mask)
    BOOST_TEST(v == 1);

  qt_locator.build(qt.root_size(), qt.nodes());

  for (uint32_t py = 0; py < root; ++py)
    for (uint32_t px = 0; px < root; ++px) {
      auto lf = qt_locator.locate(px + 0.5, py + 0.5);
      BOOST_TEST(px >= lf.x);
      BOOST_TEST(px < lf.x + lf.size);
      BOOST_TEST(py >= lf.y);
      BOOST_TEST(py < lf.y + lf.size);
    }
}

BOOST_AUTO_TEST_SUITE_END()
