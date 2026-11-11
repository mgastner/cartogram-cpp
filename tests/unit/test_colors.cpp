#define BOOST_TEST_MODULE test_colors
#include "colors.hpp"
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(ColorTests)

const Color white(255, 255, 255);

BOOST_AUTO_TEST_CASE(HexCode)
{
  BOOST_CHECK(Color("#e74c3c") == Color(0xe7, 0x4c, 0x3c));
  BOOST_CHECK(Color("#E74C3C") == Color(0xe7, 0x4c, 0x3c));
}

BOOST_AUTO_TEST_CASE(ShortHexCode)
{
  // "#abc" stands for "#aabbcc", as in CSS
  BOOST_CHECK(Color("#abc") == Color(0xaa, 0xbb, 0xcc));
}

BOOST_AUTO_TEST_CASE(HtmlColorName)
{
  BOOST_CHECK(Color("CornflowerBlue") == Color(100, 149, 237));
}

BOOST_AUTO_TEST_CASE(RgbValues)
{
  BOOST_CHECK(Color("rgb(10, 20, 30)") == Color(10, 20, 30));
}

BOOST_AUTO_TEST_CASE(InvalidHexCodesBecomeWhite)
{
  // Before, std::stoi() threw for the first two strings, which aborted the
  // program
  BOOST_CHECK(Color("#zzzzzz") == white);
  BOOST_CHECK(Color("#ffffffffffff") == white);
  BOOST_CHECK(Color("#12") == white);
  BOOST_CHECK(Color("#") == white);
}

BOOST_AUTO_TEST_SUITE_END()
