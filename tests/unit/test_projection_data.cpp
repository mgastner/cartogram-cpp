#define BOOST_TEST_MODULE test_projection_data
#include "projection_data.hpp"
#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include <cstdint>
#include <random>
#include <unordered_set>
#include <vector>

namespace
{
inline uint64_t pack_xy(uint32_t x, uint32_t y)
{
  return (uint64_t(x) << 32) | uint64_t(y);
}

std::vector<QuadtreeCorner> make_unique_random_keys(
  uint32_t lx,
  uint32_t ly,
  size_t n,
  std::mt19937 &rng)
{
  BOOST_REQUIRE(lx > 0 && ly > 0);
  const uint64_t cap = uint64_t(lx) * uint64_t(ly);
  BOOST_REQUIRE_MESSAGE(n <= cap, "n exceeds grid capacity");

  std::uniform_int_distribution<uint32_t> dx(0, lx - 1);
  std::uniform_int_distribution<uint32_t> dy(0, ly - 1);

  std::unordered_set<uint64_t> seen;
  std::vector<QuadtreeCorner> out;
  out.reserve(n);
  while (out.size() < n) {
    uint32_t x = dx(rng), y = dy(rng);
    uint64_t key = pack_xy(x, y);
    if (seen.insert(key).second)
      out.emplace_back(x, y);
  }
  return out;
}

void fill_projection(std::vector<Point> &proj, size_t n)
{
  proj.resize(n);
  for (size_t i = 0; i < n; ++i) {

    // Any deterministic, unique values
    proj[i] = Point(double(i) + 0.25, double(i) + 0.75);
  }
}
};  // namespace

BOOST_AUTO_TEST_CASE(zero_initialized_buffer_makes_unset_cells_invalid)
{
  ProjectionData pd;
  pd.reserve(8, 8);
  for (uint32_t x = 0; x < 8; ++x) {
    for (uint32_t y = 0; y < 8; ++y) {
      BOOST_TEST(!pd.is_valid_corner(x, y));
    }
  }
}

BOOST_AUTO_TEST_CASE(build_and_query_basic_mapping_and_get)
{
  const uint32_t lx = 16, ly = 12;
  ProjectionData pd;
  pd.reserve(lx, ly);

  std::vector<QuadtreeCorner> keys{{0, 0}, {5, 7}, {15, 11}, {1, 2}, {9, 3}};
  fill_projection(pd.get_projection(), keys.size());

  pd.build_fast_indexing(keys);

  for (uint32_t i = 0; i < keys.size(); ++i) {
    auto [x, y] = std::pair<uint32_t, uint32_t>(keys[i].x(), keys[i].y());
    BOOST_TEST(pd.is_valid_corner(x, y));
    BOOST_TEST(pd.offset(x, y) == i);
    BOOST_TEST(pd.get(x, y) == pd.get_projection()[i]);
  }

  // Unset cells should remain invalid
  BOOST_TEST(!pd.is_valid_corner(10, 10));
  BOOST_TEST(!pd.is_valid_corner(20, 20));
  BOOST_TEST(!pd.is_valid_corner(0, 12));
  BOOST_TEST(!pd.is_valid_corner(16, 0));
  BOOST_TEST(!pd.is_valid_corner(15, 12));
  BOOST_TEST(!pd.is_valid_corner(16, 12));
  BOOST_TEST(!pd.is_valid_corner(0, 11));
}

BOOST_AUTO_TEST_CASE(duplicate_keys_last_wins)
{
  const uint32_t lx = 8, ly = 8;
  ProjectionData pd;
  pd.reserve(lx, ly);

  // Same corner appears twice; last occurrence should win (index = 2)
  std::vector<QuadtreeCorner> keys{{1, 1}, {2, 2}, {1, 1}};
  fill_projection(pd.get_projection(), keys.size());
  pd.build_fast_indexing(keys);

  BOOST_TEST(pd.is_valid_corner(1, 1));
  BOOST_TEST(pd.offset(1, 1) == 2u);

  BOOST_TEST(pd.is_valid_corner(2, 2));
  BOOST_TEST(pd.offset(2, 2) == 1u);
}

BOOST_AUTO_TEST_CASE(generation_invalidation_and_rebuild)
{
  const uint32_t lx = 10, ly = 10;
  ProjectionData pd;
  pd.reserve(lx, ly);

  std::vector<QuadtreeCorner> keys1{{0, 0}, {1, 1}, {2, 2}};
  std::vector<QuadtreeCorner> keys2{{3, 3}, {4, 4}, {5, 5}};

  fill_projection(pd.get_projection(), keys1.size());
  pd.build_fast_indexing(keys1);

  for (auto &k : keys1) {
    BOOST_TEST(pd.is_valid_corner(k.x(), k.y()));
  }
  // keys2 should not be valid yet
  for (auto &k : keys2) {
    BOOST_TEST(!pd.is_valid_corner(k.x(), k.y()));
  }

  // Rebuild with a disjoint set; prior cells must be invalid now
  fill_projection(pd.get_projection(), keys2.size());
  pd.build_fast_indexing(keys2);

  for (auto &k : keys1) {
    BOOST_TEST(!pd.is_valid_corner(k.x(), k.y()));
  }

  for (uint32_t i = 0; i < keys2.size(); ++i) {
    BOOST_TEST(pd.is_valid_corner(keys2[i].x(), keys2[i].y()));
    BOOST_TEST(pd.offset(keys2[i].x(), keys2[i].y()) == i);
    BOOST_TEST(pd.get(keys2[i].x(), keys2[i].y()) == pd.get_projection()[i]);
  }
}

BOOST_AUTO_TEST_CASE(reserve_updates_shape_even_without_realloc)
{
  // This specifically catches any old bug where lx_/ly_ weren't updated if
  // capacity sufficed. Strategy:
  // 1) Reserve larger, build at (1,1) so slot at index 1*old_ly+1 is written.
  // 2) Reserve to smaller ly but still within existing capacity -> no
  // reallocation.
  // 3) Without rebuilding, is_valid_corner(1,1) must be false because the new
  // stride maps (1,1) to a different slot that should still be zero.

  const uint32_t old_lx = 3, old_ly = 10;  // area 30
  const uint32_t new_lx = 3, new_ly = 7;  // area 21 (<= capacity)

  ProjectionData pd;
  pd.reserve(old_lx, old_ly);

  std::vector<QuadtreeCorner> keys{{1, 1}};
  fill_projection(pd.get_projection(), keys.size());
  pd.build_fast_indexing(keys);
  BOOST_TEST(pd.is_valid_corner(1, 1));  // sanity

  // Change shape only; capacity stays
  pd.reserve(new_lx, new_ly);

  // If ly_ wasn't updated, this would still read the previously written slot
  // and be true.
  BOOST_TEST(!pd.is_valid_corner(1, 1));
}

BOOST_AUTO_TEST_CASE(copy_constructor_and_assignment_are_deep)
{
  const uint32_t lx = 12, ly = 9;
  ProjectionData a;
  a.reserve(lx, ly);

  std::vector<QuadtreeCorner> keys{{2, 3}, {4, 5}, {6, 7}};
  fill_projection(a.get_projection(), keys.size());
  a.build_fast_indexing(keys);

  // Copy construct
  ProjectionData b(a);

  // Both should see same validity on their own data
  for (uint32_t i = 0; i < keys.size(); ++i) {
    uint32_t x = keys[i].x(), y = keys[i].y();
    BOOST_TEST(a.is_valid_corner(x, y));
    BOOST_TEST(b.is_valid_corner(x, y));
    BOOST_TEST(a.get(x, y) == b.get(x, y));
  }

  // Mutate 'a' with a new generation; 'b' must remain unaffected.
  std::vector<QuadtreeCorner> keys2{{0, 1}, {1, 2}};
  fill_projection(a.get_projection(), keys2.size());
  a.build_fast_indexing(keys2);

  for (auto &k : keys) {
    BOOST_TEST(!a.is_valid_corner(k.x(), k.y()));  // invalidated in 'a'
    BOOST_TEST(b.is_valid_corner(k.x(), k.y()));  // still valid in 'b'
  }

  // Copy assign from 'a' into a fresh object and verify it matches 'a'
  ProjectionData c;
  c = a;
  for (auto &k : keys2) {
    BOOST_TEST(c.is_valid_corner(k.x(), k.y()));
    BOOST_TEST(c.get(k.x(), k.y()) == a.get(k.x(), k.y()));
  }
}

BOOST_AUTO_TEST_CASE(move_constructor_and_assignment_transfer_state)
{
  const uint32_t lx = 20, ly = 15;
  ProjectionData src;
  src.reserve(lx, ly);

  std::vector<QuadtreeCorner> keys{{7, 3}, {10, 14}, {19, 0}};
  fill_projection(src.get_projection(), keys.size());
  src.build_fast_indexing(keys);

  // Move construct
  ProjectionData dst(std::move(src));
  for (uint32_t i = 0; i < keys.size(); ++i) {
    uint32_t x = keys[i].x(), y = keys[i].y();
    BOOST_TEST(dst.is_valid_corner(x, y));
    BOOST_TEST(dst.offset(x, y) == i);
    BOOST_TEST(dst.get(x, y) == dst.get_projection()[i]);
  }

  // Move assign
  ProjectionData dst2;
  dst2 = std::move(dst);
  for (uint32_t i = 0; i < keys.size(); ++i) {
    uint32_t x = keys[i].x(), y = keys[i].y();
    BOOST_TEST(dst2.is_valid_corner(x, y));
    BOOST_TEST(dst2.offset(x, y) == i);
    BOOST_TEST(dst2.get(x, y) == dst2.get_projection()[i]);
  }
}

BOOST_AUTO_TEST_CASE(num_unique_corners_matches_projection_size)
{
  ProjectionData pd;
  pd.reserve(5, 5);

  for (size_t n = 0; n < 20; ++n) {
    pd.get_projection().resize(n, Point(0.0, 0.0));
    BOOST_TEST(pd.num_unique_corners() == n);
  }
}

BOOST_AUTO_TEST_CASE(quadtreecorner_ops_and_conversion)
{
  QuadtreeCorner a(1, 2), b(1, 2), c(2, 1);
  BOOST_TEST(a == b);
  BOOST_TEST(!(a == c));

  // Implicit conversion to Point should compile and yield the same as
  // Point(x,y)
  Point p1 = static_cast<Point>(a);
  Point p2(1.0, 2.0);
  BOOST_TEST(p1 == p2);
}

BOOST_AUTO_TEST_CASE(randomized_fuzz_roundtrips)
{
  const uint32_t lx = 64, ly = 96;
  ProjectionData pd;
  pd.reserve(lx, ly);

  std::mt19937 rng(123456789);

  // Try several rounds with different sizes
  for (std::size_t round = 0; round < 5; ++round) {
    constexpr std::size_t BASE = 500;
    constexpr std::size_t STEP = 100;
    const std::size_t n = BASE + round * STEP;
    auto keys = make_unique_random_keys(lx, ly, n, rng);

    fill_projection(pd.get_projection(), n);
    pd.build_fast_indexing(keys);

    // Sample a subset and verify roundtrip (is_valid -> offset -> get)
    std::uniform_int_distribution<size_t> pick(0, n - 1);
    for (int t = 0; t < 500; ++t) {
      const auto &k = keys[pick(rng)];
      BOOST_TEST(pd.is_valid_corner(k.x(), k.y()));
      uint32_t idx = pd.offset(k.x(), k.y());
      BOOST_REQUIRE(idx < pd.get_projection().size());
      BOOST_TEST(pd.get(k.x(), k.y()) == pd.get_projection()[idx]);
    }

    // Pick random cells *not* in the set and ensure they are invalid
    std::unordered_set<uint64_t> S;
    S.reserve(keys.size());
    for (auto &k : keys)
      S.insert(pack_xy(k.x(), k.y()));

    std::uniform_int_distribution<uint32_t> dx(0, lx - 1), dy(0, ly - 1);
    int checked = 0;
    while (checked < 200) {
      uint32_t x = dx(rng), y = dy(rng);
      if (S.count(pack_xy(x, y)) == 0) {
        BOOST_TEST(!pd.is_valid_corner(x, y));
        ++checked;
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(capacity_growth_zero_initializes_new_regions)
{
  // Grow capacity; new slots must be invalid (zeroed)
  ProjectionData pd;
  pd.reserve(8, 8);

  std::vector<QuadtreeCorner> keys{{2, 2}, {3, 3}};
  fill_projection(pd.get_projection(), keys.size());
  pd.build_fast_indexing(keys);

  // Grow to a larger shape (forces reallocation/zero-init)
  pd.reserve(32, 16);

  // Previously unused cell in the newly reachable area must be invalid
  BOOST_TEST(!pd.is_valid_corner(31, 15));
  BOOST_TEST(!pd.is_valid_corner(0, 15));
  BOOST_TEST(!pd.is_valid_corner(31, 0));
}
