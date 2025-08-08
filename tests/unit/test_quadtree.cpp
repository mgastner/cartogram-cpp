#define BOOST_TEST_MODULE test_quadtree
#include "quadtree.hpp"
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

BOOST_AUTO_TEST_SUITE_END()
