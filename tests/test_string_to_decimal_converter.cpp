#define BOOST_TEST_MODULE StringToDecimalConverterTest
#include "string_to_decimal_converter.hpp"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(TestIsStrNA_True)
{
  BOOST_CHECK(StringToDecimalConverter::is_str_NA("NA"));
}

BOOST_AUTO_TEST_CASE(TestIsStrNA_False)
{
  BOOST_CHECK(!StringToDecimalConverter::is_str_NA("123"));
  BOOST_CHECK(!StringToDecimalConverter::is_str_NA(""));
}

BOOST_AUTO_TEST_CASE(TestValidCharacters_BasicNumeric)
{
  BOOST_CHECK(StringToDecimalConverter::is_str_valid_characters("123456789"));
}

BOOST_AUTO_TEST_CASE(TestValidCharacters_CommasAndPoints)
{
  BOOST_CHECK(StringToDecimalConverter::is_str_valid_characters("123,456.789"));
}

BOOST_AUTO_TEST_CASE(TestValidCharacters_InvalidCharacters)
{
  BOOST_CHECK(!StringToDecimalConverter::is_str_valid_characters("abc"));
}

BOOST_AUTO_TEST_CASE(TestValidCharacters_NegativeNumber)
{
  BOOST_CHECK(StringToDecimalConverter::is_str_valid_characters("-123.456"));
}

BOOST_AUTO_TEST_CASE(TestValidCharacters_MultipleHyphens)
{
  BOOST_CHECK(!StringToDecimalConverter::is_str_valid_characters("--123.456"));
}

BOOST_AUTO_TEST_CASE(TestValidCharacters_HyphenNotAtStart)
{
  BOOST_CHECK(!StringToDecimalConverter::is_str_valid_characters("123-456.789"));
}

BOOST_AUTO_TEST_CASE(TestCorrectFormat_SimpleNumeric)
{
  BOOST_CHECK(StringToDecimalConverter::is_str_correct_format("123456789"));
}

BOOST_AUTO_TEST_CASE(TestCorrectFormat_CommaAndPoint)
{
  BOOST_CHECK(StringToDecimalConverter::is_str_correct_format("123,456.789"));
}

BOOST_AUTO_TEST_CASE(TestCorrectFormat_PointAndComma)
{
  BOOST_CHECK(StringToDecimalConverter::is_str_correct_format("123.456,789"));
}

BOOST_AUTO_TEST_CASE(TestCorrectFormat_MultipleCommasAndPoints)
{
  BOOST_CHECK(
    !StringToDecimalConverter::is_str_correct_format("123,456.789,123"));
}

BOOST_AUTO_TEST_CASE(TestCorrectFormat_CommaAtEnd)
{
  BOOST_CHECK(!StringToDecimalConverter::is_str_correct_format("123456789,"));
}

BOOST_AUTO_TEST_CASE(TestCorrectFormat_PointAtEnd)
{
  BOOST_CHECK(!StringToDecimalConverter::is_str_correct_format("123456789."));
}

BOOST_AUTO_TEST_CASE(TestParseStrWithCommasAndPoints)
{
  // Assuming commas are thousands separators and points are decimal separators
  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("123,456.789"),
    123456.789);

  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("123,123,456.789"),
    123123456.789);

  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("123,123,123,456.789"),
    123123123456.789);

  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("12,31,23,12,34,56.789"),
    123123123456.789);

  // Assuming points are thousands separators and commas are decimal separators
  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("123.456,789"),
    123456.789);

  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("123.123.456,789"),
    123123456.789);

  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("123.123.123.456,789"),
    123123123456.789);

  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("12.31.23.12.34.56,789"),
    123123123456.789);
}

BOOST_AUTO_TEST_CASE(TestParseStrWithOnlyCommas)
{
  // Single comma acting as a decimal point
  BOOST_CHECK_EQUAL(StringToDecimalConverter::parse_str("123456,78"), 123456.78);

  // Comma as thousands separator
  BOOST_CHECK_EQUAL(StringToDecimalConverter::parse_str("123,456"), 123456);

  // Multiple commas as thousands separators
  BOOST_CHECK_EQUAL(StringToDecimalConverter::parse_str("1,234,567"), 1234567);

  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("123,456,789,123"),
    123456789123);
}

BOOST_AUTO_TEST_CASE(TestParseStrWithOnlyPoints)
{
  // Single point as decimal separator, otherwise there would be another point
  // after 1
  BOOST_CHECK_EQUAL(StringToDecimalConverter::parse_str("1234.567"), 1234.567);

  // Single point as thousands separator
  BOOST_CHECK_EQUAL(StringToDecimalConverter::parse_str("12.345"), 12345);

  // Multiple points as thousands separators
  BOOST_CHECK_EQUAL(StringToDecimalConverter::parse_str("1.234.567"), 1234567);

  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("123.456.789.123"),
    123456789123);
}

BOOST_AUTO_TEST_CASE(TestParseStrNegativeValue)
{
  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("-123456789"),
    -123456789);
}
