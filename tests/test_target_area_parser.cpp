#define BOOST_TEST_MODULE TargetAreaParserTest
#include "target_area_parser.h"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(TestIsAreaStrNA_True)
{
  BOOST_CHECK(TargetAreaParser::is_area_str_NA("NA"));
}

BOOST_AUTO_TEST_CASE(TestIsAreaStrNA_False)
{
  BOOST_CHECK(!TargetAreaParser::is_area_str_NA("123"));
  BOOST_CHECK(!TargetAreaParser::is_area_str_NA(""));
}

BOOST_AUTO_TEST_CASE(TestValidCharacters_BasicNumeric)
{
  BOOST_CHECK(TargetAreaParser::is_area_str_valid_characters("123456789"));
}

BOOST_AUTO_TEST_CASE(TestValidCharacters_CommasAndPoints)
{
  BOOST_CHECK(TargetAreaParser::is_area_str_valid_characters("123,456.789"));
}

BOOST_AUTO_TEST_CASE(TestValidCharacters_InvalidCharacters)
{
  BOOST_CHECK(!TargetAreaParser::is_area_str_valid_characters("abc"));
}

BOOST_AUTO_TEST_CASE(TestValidCharacters_NegativeNumber)
{
  BOOST_CHECK(TargetAreaParser::is_area_str_valid_characters("-123.456"));
}

BOOST_AUTO_TEST_CASE(TestValidCharacters_MultipleHyphens)
{
  BOOST_CHECK(!TargetAreaParser::is_area_str_valid_characters("--123.456"));
}

BOOST_AUTO_TEST_CASE(TestValidCharacters_HyphenNotAtStart)
{
  BOOST_CHECK(!TargetAreaParser::is_area_str_valid_characters("123-456.789"));
}

BOOST_AUTO_TEST_CASE(TestCorrectFormat_SimpleNumeric)
{
  BOOST_CHECK(TargetAreaParser::is_area_str_correct_format("123456789"));
}

BOOST_AUTO_TEST_CASE(TestCorrectFormat_CommaAndPoint)
{
  BOOST_CHECK(TargetAreaParser::is_area_str_correct_format("123,456.789"));
}

BOOST_AUTO_TEST_CASE(TestCorrectFormat_PointAndComma)
{
  BOOST_CHECK(TargetAreaParser::is_area_str_correct_format("123.456,789"));
}

BOOST_AUTO_TEST_CASE(TestCorrectFormat_MultipleCommasAndPoints)
{
  BOOST_CHECK(
    !TargetAreaParser::is_area_str_correct_format("123,456.789,123"));
}

BOOST_AUTO_TEST_CASE(TestCorrectFormat_CommaAtEnd)
{
  BOOST_CHECK(!TargetAreaParser::is_area_str_correct_format("123456789,"));
}

BOOST_AUTO_TEST_CASE(TestCorrectFormat_PointAtEnd)
{
  BOOST_CHECK(!TargetAreaParser::is_area_str_correct_format("123456789."));
}

BOOST_AUTO_TEST_CASE(TestParseAreaStrWithCommasAndPoints)
{
  // Assuming commas are thousands separators and points are decimal separators
  BOOST_CHECK_EQUAL(
    TargetAreaParser::parse_area_str("123,456.789"),
    123456.789);

  BOOST_CHECK_EQUAL(
    TargetAreaParser::parse_area_str("123,123,456.789"),
    123123456.789);

  BOOST_CHECK_EQUAL(
    TargetAreaParser::parse_area_str("123,123,123,456.789"),
    123123123456.789);

  BOOST_CHECK_EQUAL(
    TargetAreaParser::parse_area_str("12,31,23,12,34,56.789"),
    123123123456.789);

  // Assuming points are thousands separators and commas are decimal separators
  BOOST_CHECK_EQUAL(
    TargetAreaParser::parse_area_str("123.456,789"),
    123456.789);

  BOOST_CHECK_EQUAL(
    TargetAreaParser::parse_area_str("123.123.456,789"),
    123123456.789);

  BOOST_CHECK_EQUAL(
    TargetAreaParser::parse_area_str("123.123.123.456,789"),
    123123123456.789);

  BOOST_CHECK_EQUAL(
    TargetAreaParser::parse_area_str("12.31.23.12.34.56,789"),
    123123123456.789);
}

BOOST_AUTO_TEST_CASE(TestParseAreaStrWithOnlyCommas)
{
  // Single comma acting as a decimal point
  BOOST_CHECK_EQUAL(TargetAreaParser::parse_area_str("123456,78"), 123456.78);

  // Comma as thousands separator
  BOOST_CHECK_EQUAL(TargetAreaParser::parse_area_str("123,456"), 123456);

  // Multiple commas as thousands separators
  BOOST_CHECK_EQUAL(TargetAreaParser::parse_area_str("1,234,567"), 1234567);

  BOOST_CHECK_EQUAL(
    TargetAreaParser::parse_area_str("123,456,789,123"),
    123456789123);
}

BOOST_AUTO_TEST_CASE(TestParseAreaStrWithOnlyPoints)
{
  // Single point as decimal separator, otherwise there would be another point
  // after 1
  BOOST_CHECK_EQUAL(TargetAreaParser::parse_area_str("1234.567"), 1234.567);

  // Single point as thousands separator
  BOOST_CHECK_EQUAL(TargetAreaParser::parse_area_str("12.345"), 12345);

  // Multiple points as thousands separators
  BOOST_CHECK_EQUAL(TargetAreaParser::parse_area_str("1.234.567"), 1234567);

  BOOST_CHECK_EQUAL(
    TargetAreaParser::parse_area_str("123.456.789.123"),
    123456789123);
}

BOOST_AUTO_TEST_CASE(TestParseAreaStrNegativeValue)
{
  BOOST_CHECK_EQUAL(
    TargetAreaParser::parse_area_str("-123456789"),
    -123456789);
}
