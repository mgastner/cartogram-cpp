/*
String to Decimal Parsing Algorithm:
*) Check validity of string: only 0-9, ',' and '.' are allowed.

*) Contains both points and commas (i.e: 123,456,789.123 or
123.456.789,123):
  i) Number of commas and points both are more than 1
  (123,456,789.123.456): This format does not belong to any known
  convention, so exit with error.
  ii) If number of one type of separator is 1 and the other is more than 1,
  we assume the separator with count 1 is the decimal separator and it must
  come before the other separator. (valid: 123,456,789.123 -> 123456789.123,
  123.456.789,123 -> 123456789.123; invalid: 123,456.789.123 -> error)
  iii) Comma appears first (123,456,789.123):
  We assume commas are used as big-number separators and remove them.
  (123,456,789.123 -> 123456789.123).
  iv) Point appears first (123.456.789,123): We assume points are used as
  big-number separators and remove them (123.456.789,123 -> 123456789.123)

*) Contains only commas (i.e: 123,456,789 or 123456,78):
  i) Only one comma is present and there are two digits after comma. We
  treat the comma as decimal separator. (123456,78 -> 123456.78)
  ii) All other cases, we assume the number does not contain any fractional
  part. So, we remove the commas and parse as usual. (123,456,789 ->
  123456789)

*) Contains only points (i.e., "123.456.789" or "123456.789"):
  i) Only one point is present:
    - If the point is followed by exactly three digits and the total length
    of the number (excluding the point) is more than four, assume the point is
a thousands separator and remove it. (e.g., "1.234" -> "1234", "12.345" ->
    "12345").
    - In other cases, assume the point is a decimal separator and keep it.
      (e.g., "123.45" remains "123.45", "1234.5" remains "1234.5").
  ii) Multiple points are present:
    - Assume all points are used as thousands separators and remove them.
      (e.g., "1.234.567" -> "1234567").
*/

#include "string_to_decimal_converter.h"
#include <algorithm>
#include <cassert>
#include <cctype>

const std::string StringToDecimalConverter::NA_ = "NA";

bool StringToDecimalConverter::is_valid_char(char ch)
{
  return (std::isdigit(ch)) || ch == point_ || ch == comma_ || ch == minus_;
}

std::string StringToDecimalConverter::remove_char(std::string str, char ch)
{
  str.erase(std::remove(str.begin(), str.end(), ch), str.end());
  return str;
}

int StringToDecimalConverter::count_char(const std::string &str, char ch)
{
  return std::count(str.begin(), str.end(), ch);
}

bool StringToDecimalConverter::has_multiple_commas_and_points(
  const std::string &str)
{
  return count_char(str, comma_) > 1 && count_char(str, point_) > 1;
}

bool StringToDecimalConverter::has_separator_at_the_end(const std::string &str)
{
  assert(str.length() > 0);
  return (str.back() == comma_ || str.back() == point_);
}

bool StringToDecimalConverter::has_invalid_comma_point_sequence(
  const std::string &str)
{
  assert(!has_multiple_commas_and_points(str));

  const int comma_count = count_char(str, comma_);
  const int point_count = count_char(str, point_);

  // If both comma and point are not found, the sequence is valid
  if (comma_count == 0 || point_count == 0) {
    return false;
  }

  if (comma_count > 1) {
    assert(point_count == 1);
    size_t last_comma_pos = str.rfind(comma_);
    size_t point_pos = str.find(point_);

    // Check if the point is the rightmost among all separators
    if (point_pos < last_comma_pos) {
      return true;
    }
  }

  if (point_count > 1) {
    assert(comma_count == 1);
    size_t last_point_pos = str.rfind(point_);
    size_t comma_pos = str.find(comma_);

    // Check if the comma is the rightmost among all separators
    if (comma_pos < last_point_pos) {
      return true;
    }
  }

  return false;
}

bool StringToDecimalConverter::is_str_NA(const std::string &str)
{
  return (str.compare(NA_) == 0);
}

bool StringToDecimalConverter::is_str_valid_characters(const std::string &str)
{
  if (str.empty()) {
    return false;
  }

  // Allow str being "NA"
  if (str == NA_) {
    return true;
  }

  // Only 0 to 9, '.', '-', and ',' are allowed
  for (const auto &c : str) {
    if (!is_valid_char(c)) {
      return false;
    }
  }

  // '-' can only be used once
  if (count_char(str, minus_) > 1) {
    return false;
  }

  // '-' can only be used at the beginning
  if (count_char(str, minus_) == 1 and str[0] != minus_) {
    return false;
  }
  return true;
}

bool StringToDecimalConverter::is_str_correct_format(const std::string &str)
{
  assert(is_str_valid_characters(str));
  assert(is_str_NA(str) == false);

  // if the number of commas and points both are more than 1, then this format
  // does not belong to any known convention
  if (has_multiple_commas_and_points(str)) {
    return false;
  }

  // Check for commas before and after a point, or points before and after a
  // comma
  if (has_invalid_comma_point_sequence(str)) {
    return false;
  }

  // Check for separators at the end of the string
  if (has_separator_at_the_end(str)) {
    return false;
  }

  return true;
}

double StringToDecimalConverter::parse_str(const std::string &str)
{
  assert(is_str_correct_format(str));
  assert(!is_str_NA(str));

  std::string processed_str = str;

  int comma_count = count_char(str, comma_);
  int point_count = count_char(str, point_);

  // Contains both commas and points
  if (comma_count > 0 && point_count > 0) {
    if (str.find(point_) > str.find(comma_)) {

      // Commas as thousand separators, remove them
      processed_str = remove_char(processed_str, comma_);
    } else {

      // Points as thousand separators, remove them
      processed_str = remove_char(processed_str, point_);

      assert(processed_str.find(comma_) != std::string::npos);

      // Replace the comma with a point
      processed_str[processed_str.find(comma_)] = point_;
    }
  }
  // Contains only commas
  else if (comma_count > 0) {

    // If only one comma and two digits after it, treat as a decimal
    if (comma_count == 1 && str.size() - str.find(comma_) == 3) {
      processed_str[processed_str.find(comma_)] = point_;
    } else {

      // Otherwise, remove all commas
      processed_str = remove_char(processed_str, comma_);
    }
  }
  // Contains only points
  else if (point_count > 0) {

    // Check for the presence of multiple points
    if (point_count > 1) {

      // If there are multiple points, assume all points are used as
      // thousands separators and remove them.
      processed_str = remove_char(processed_str, point_);
    } else {

      // If there is only one point, check the number of digits following it.
      size_t point_pos = processed_str.find(point_);
      size_t digits_after_point = processed_str.length() - point_pos - 1;

      if (
        digits_after_point == 3 && processed_str.length() > 4 &&
        processed_str.size() < 8) {

        // If exactly three digits follow the point and the total length of
        // the number (excluding the point) is more than four, and
        // less than 8, assume the point is a thousands separator and remove
        // it. (e.g., "1.234" -> "1234", "123.456" -> "123456").
        processed_str = remove_char(processed_str, point_);
      }

      // In other cases, assume the point is a decimal separator and keep it.
      // (e.g., "
    }

    // Parse the str or return it as needed.
  }

  return std::stod(processed_str);
}
