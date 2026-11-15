#define BOOST_TEST_MODULE test_blur_width
#include "inset_state.hpp"
#include <boost/test/included/unit_test.hpp>

// Before the first integration, the blur width is 1/128 of the long side of
// the grid, but at least 1 cell. Each failed integration doubles the blur
// width. Long sides of fewer than 256 cells, which only the option -n
// produces, must give a well-defined blur width as well.

namespace
{

double blur_width_of_grid(
  const unsigned int lx,
  const unsigned int ly,
  const unsigned int n_fails = 0)
{
  Arguments args{};
  InsetState inset("C", args);
  inset.set_grid_dimensions(lx, ly);
  for (unsigned int i = 0; i < n_fails; ++i) {
    inset.increment_n_fails_during_flatten_density();
  }
  return inset.blur_width();
}

}  // namespace

BOOST_AUTO_TEST_CASE(long_side_below_256_cells)
{
  BOOST_TEST(blur_width_of_grid(128, 64) == 1.0);
  BOOST_TEST(blur_width_of_grid(64, 128) == 1.0);
  BOOST_TEST(blur_width_of_grid(128, 128) == 1.0);
  BOOST_TEST(blur_width_of_grid(64, 32) == 1.0);
  BOOST_TEST(blur_width_of_grid(32, 64) == 1.0);
}

BOOST_AUTO_TEST_CASE(long_side_of_256_cells_or_more)
{
  BOOST_TEST(blur_width_of_grid(256, 128) == 2.0);
  BOOST_TEST(blur_width_of_grid(128, 256) == 2.0);
  BOOST_TEST(blur_width_of_grid(512, 512) == 4.0);
  BOOST_TEST(blur_width_of_grid(1024, 256) == 8.0);
  BOOST_TEST(blur_width_of_grid(256, 2048) == 16.0);
}

BOOST_AUTO_TEST_CASE(each_failure_doubles_the_blur_width)
{
  BOOST_TEST(blur_width_of_grid(64, 32, 1) == 2.0);
  BOOST_TEST(blur_width_of_grid(128, 64, 2) == 4.0);
  BOOST_TEST(blur_width_of_grid(256, 128, 1) == 4.0);
  BOOST_TEST(blur_width_of_grid(512, 256, 3) == 32.0);
}
