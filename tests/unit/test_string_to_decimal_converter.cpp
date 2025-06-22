#define BOOST_TEST_MODULE test_string_to_decimal_converter
#include "string_to_decimal_converter.hpp"
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(StringToDecimalConverterTests)

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
  BOOST_CHECK(
    StringToDecimalConverter::is_str_valid_characters("123,456.789"));
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
  BOOST_CHECK(
    !StringToDecimalConverter::is_str_valid_characters("123-456.789"));
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

BOOST_AUTO_TEST_CASE(TestValidCharacters_ScientificNotation)
{
  // Basic scientific notation
  BOOST_CHECK(StringToDecimalConverter::is_str_valid_characters("1.23e4"));
  BOOST_CHECK(StringToDecimalConverter::is_str_valid_characters("1.23E4"));

  // Negative exponents
  BOOST_CHECK(StringToDecimalConverter::is_str_valid_characters("1.23e-4"));
  BOOST_CHECK(StringToDecimalConverter::is_str_valid_characters("1.23E-4"));

  // With commas as thousand separators
  BOOST_CHECK(StringToDecimalConverter::is_str_valid_characters("1,234.56e4"));
  BOOST_CHECK(StringToDecimalConverter::is_str_valid_characters("1.234,56E4"));

  // Negative numbers
  BOOST_CHECK(StringToDecimalConverter::is_str_valid_characters("-1.23e4"));
  BOOST_CHECK(StringToDecimalConverter::is_str_valid_characters("-1.23E-4"));

  // Invalid scientific notation
  BOOST_CHECK(!StringToDecimalConverter::is_str_valid_characters(
    "1.23e"));  // Missing exponent
  BOOST_CHECK(!StringToDecimalConverter::is_str_valid_characters(
    "e4"));  // Missing mantissa
  BOOST_CHECK(!StringToDecimalConverter::is_str_valid_characters(
    "1.23ee4"));  // Multiple e's
  BOOST_CHECK(!StringToDecimalConverter::is_str_valid_characters(
    "1.23e4.5"));  // Non-integer exponent
  BOOST_CHECK(!StringToDecimalConverter::is_str_valid_characters(
    "1.23e-"));  // Incomplete negative exponent
}

BOOST_AUTO_TEST_CASE(TestCorrectFormat_ScientificNotation)
{
  // Valid formats
  BOOST_CHECK(StringToDecimalConverter::is_str_correct_format("1.23e4"));
  BOOST_CHECK(StringToDecimalConverter::is_str_correct_format("1.23E4"));
  BOOST_CHECK(StringToDecimalConverter::is_str_correct_format("1,234.56e4"));
  BOOST_CHECK(StringToDecimalConverter::is_str_correct_format("1.234,56E4"));
  BOOST_CHECK(StringToDecimalConverter::is_str_correct_format("-1.23e4"));
  BOOST_CHECK(StringToDecimalConverter::is_str_correct_format("-1.23E-4"));

  // Invalid formats
  BOOST_CHECK(!StringToDecimalConverter::is_str_correct_format(
    "1.23e4.5"));  // Non-integer exponent
  BOOST_CHECK(!StringToDecimalConverter::is_str_correct_format(
    "1.23e4,"));  // Comma at end
  BOOST_CHECK(!StringToDecimalConverter::is_str_correct_format(
    "1.23e4."));  // Point at end
  BOOST_CHECK(!StringToDecimalConverter::is_str_correct_format(
    "1.23,456.789e4"));  // Multiple separators
  BOOST_CHECK(!StringToDecimalConverter::is_str_correct_format(
    "1.23e4e5"));  // Multiple e's
  BOOST_CHECK(!StringToDecimalConverter::is_str_correct_format(
    "1.23e"));  // Missing exponent
  BOOST_CHECK(!StringToDecimalConverter::is_str_correct_format(
    "e4"));  // Missing mantissa
}

BOOST_AUTO_TEST_CASE(TestParseStr_ScientificNotation)
{
  // Test parsing with point as decimal separator
  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("1.23e4", true),
    "1.23e4");
  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("1,234.56e4", true),
    "1234.56e4");
  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("-1.23e-4", true),
    "-1.23e-4");

  // Test parsing with comma as decimal separator
  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("1,23e4", false),
    "1.23e4");
  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("1.234,56E4", false),
    "1234.56E4");
  BOOST_CHECK_EQUAL(
    StringToDecimalConverter::parse_str("-1,23e-4", false),
    "-1.23e-4");
}

BOOST_AUTO_TEST_SUITE_END()
