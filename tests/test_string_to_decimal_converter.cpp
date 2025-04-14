#include "string_to_decimal_converter.hpp"
#include <gtest/gtest.h>

TEST(StringToDecimalConverterTest, TestIsStrNA_True)
{
  EXPECT_TRUE(StringToDecimalConverter::is_str_NA("NA"));
}

TEST(StringToDecimalConverterTest, TestIsStrNA_False)
{
  EXPECT_FALSE(StringToDecimalConverter::is_str_NA("123"));
  EXPECT_FALSE(StringToDecimalConverter::is_str_NA(""));
}

TEST(StringToDecimalConverterTest, TestValidCharacters_BasicNumeric)
{
  EXPECT_TRUE(StringToDecimalConverter::is_str_valid_characters("123456789"));
}

TEST(StringToDecimalConverterTest, TestValidCharacters_CommasAndPoints)
{
  EXPECT_TRUE(
    StringToDecimalConverter::is_str_valid_characters("123,456.789"));
}

TEST(StringToDecimalConverterTest, TestValidCharacters_InvalidCharacters)
{
  EXPECT_FALSE(StringToDecimalConverter::is_str_valid_characters("abc"));
}

TEST(StringToDecimalConverterTest, TestValidCharacters_NegativeNumber)
{
  EXPECT_TRUE(StringToDecimalConverter::is_str_valid_characters("-123.456"));
}

TEST(StringToDecimalConverterTest, TestValidCharacters_MultipleHyphens)
{
  EXPECT_FALSE(StringToDecimalConverter::is_str_valid_characters("--123.456"));
}

TEST(StringToDecimalConverterTest, TestValidCharacters_HyphenNotAtStart)
{
  EXPECT_FALSE(
    StringToDecimalConverter::is_str_valid_characters("123-456.789"));
}

TEST(StringToDecimalConverterTest, TestCorrectFormat_SimpleNumeric)
{
  EXPECT_TRUE(StringToDecimalConverter::is_str_correct_format("123456789"));
}

TEST(StringToDecimalConverterTest, TestCorrectFormat_CommaAndPoint)
{
  EXPECT_TRUE(StringToDecimalConverter::is_str_correct_format("123,456.789"));
}

TEST(StringToDecimalConverterTest, TestCorrectFormat_PointAndComma)
{
  EXPECT_TRUE(StringToDecimalConverter::is_str_correct_format("123.456,789"));
}

TEST(StringToDecimalConverterTest, TestCorrectFormat_MultipleCommasAndPoints)
{
  EXPECT_FALSE(
    StringToDecimalConverter::is_str_correct_format("123,456.789,123"));
}

TEST(StringToDecimalConverterTest, TestCorrectFormat_CommaAtEnd)
{
  EXPECT_FALSE(StringToDecimalConverter::is_str_correct_format("123456789,"));
}

TEST(StringToDecimalConverterTest, TestCorrectFormat_PointAtEnd)
{
  EXPECT_FALSE(StringToDecimalConverter::is_str_correct_format("123456789."));
}
