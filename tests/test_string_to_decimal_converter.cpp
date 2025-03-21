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
