#include "string_to_decimal_converter.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>

const std::string StringToDecimalConverter::NA_ = "NA";

bool StringToDecimalConverter::is_valid_char(char ch)
{
  return (std::isdigit(ch)) || ch == point_ || ch == comma_ || ch == minus_ || ch == 'e' || ch == 'E';
}

std::string StringToDecimalConverter::remove_char(std::string str, char ch)
{
  str.erase(std::remove(str.begin(), str.end(), ch), str.end());
  return str;
}

static unsigned int count_char(const std::string &str, char ch)
{
  return static_cast<unsigned int>(std::count(str.begin(), str.end(), ch));
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

  const unsigned int comma_count = count_char(str, comma_);
  const unsigned int point_count = count_char(str, point_);

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

  // Only 0 to 9, '.', '-', ',', 'e', and 'E' are allowed
  for (const auto &c : str) {
    if (!is_valid_char(c)) {
      return false;
    }
  }

  // '-' can only be used once in the mantissa and once in the exponent
  size_t first_e = str.find_first_of("eE");
  if (first_e != std::string::npos) {
    // Check mantissa part
    std::string mantissa = str.substr(0, first_e);
    if (count_char(mantissa, minus_) > 1) {
      return false;
    }
    if (count_char(mantissa, minus_) == 1 && mantissa[0] != minus_) {
      return false;
    }
    // Check exponent part
    std::string exponent = str.substr(first_e + 1);
    if (count_char(exponent, minus_) > 1) {
      return false;
    }
    if (count_char(exponent, minus_) == 1 && exponent[0] != minus_) {
      return false;
    }
  } else {
    // No scientific notation, check as before
    if (count_char(str, minus_) > 1) {
      return false;
    }
    if (count_char(str, minus_) == 1 && str[0] != minus_) {
      return false;
    }
  }

  // Check for valid scientific notation format
  size_t e_count = count_char(str, 'e') + count_char(str, 'E');
  if (e_count > 1) {
    return false;
  }
  if (e_count == 1) {
    size_t e_pos = str.find_first_of("eE");
    // Must have digits before and after 'e'/'E'
    if (e_pos == 0 || e_pos == str.length() - 1) {
      return false;
    }
    // Must have at least one digit after 'e'/'E'
    bool has_digit_after = false;
    for (size_t i = e_pos + 1; i < str.length(); i++) {
      if (std::isdigit(str[i])) {
        has_digit_after = true;
        break;
      }
    }
    if (!has_digit_after) {
      return false;
    }
  }

  return true;
}

bool StringToDecimalConverter::is_str_correct_format(const std::string &str)
{
  assert(is_str_valid_characters(str));
  assert(is_str_NA(str) == false);

  // Handle scientific notation separately
  size_t e_pos = str.find_first_of("eE");
  if (e_pos != std::string::npos) {
    // Check mantissa part
    std::string mantissa = str.substr(0, e_pos);
    if (has_multiple_commas_and_points(mantissa)) {
      return false;
    }
    if (has_invalid_comma_point_sequence(mantissa)) {
      return false;
    }
    if (has_separator_at_the_end(mantissa)) {
      return false;
    }
    // Exponent part should only contain digits and optional minus sign
    std::string exponent = str.substr(e_pos + 1);
    for (char c : exponent) {
      if (!std::isdigit(c) && c != minus_) {
        return false;
      }
    }
    return true;
  }

  // Original format validation for non-scientific notation
  if (has_multiple_commas_and_points(str)) {
    return false;
  }
  if (has_invalid_comma_point_sequence(str)) {
    return false;
  }
  if (has_separator_at_the_end(str)) {
    return false;
  }

  return true;
}

static bool is_comma_as_decimal_separator(const std::string &str)
{
  unsigned int comma_count = count_char(str, ',');
  unsigned int point_count = count_char(str, '.');

  // Case 1: One comma, and if points exist, the last point appears before the
  // comma
  if (comma_count == 1) {
    size_t comma_pos = str.find(',');
    size_t last_point_pos = str.rfind('.');

    // Check if the last point (if any) appears before the comma
    if (point_count > 0 && last_point_pos < comma_pos) {
      return true;
    }

    // Case 2: One comma, no point, and digits after the comma differ from 3
    if (point_count == 0) {
      size_t digits_after_comma = str.size() - comma_pos - 1;
      if (digits_after_comma != 3) {
        return true;
      }
    }
  }

  // Case 3: If comma count is 0 and there are more than 1 point, assume
  // comma as a decimal separator
  if (comma_count == 0 && point_count > 1) {
    return true;
  }

  return false;
}

static bool is_point_as_decimal_separator(const std::string &str)
{
  unsigned int comma_count = count_char(str, ',');
  unsigned int point_count = count_char(str, '.');

  // Case 1: One point, and if commas exist, the last comma appears before the
  // point
  if (point_count == 1) {
    size_t point_pos = str.find('.');
    size_t last_comma_pos = str.rfind(',');

    // Check if the last comma (if any) appears before the point
    if (comma_count > 0 && last_comma_pos < point_pos) {
      return true;
    }

    // Case 2: One point, no comma, and digits after the point differ from 3
    if (comma_count == 0) {
      size_t digits_after_point = str.size() - point_pos - 1;

      if (digits_after_point != 3) {
        return true;
      }
    }
  }

  // Case 3: If point count is 0 and there are more than 1 comma, assume point
  // as a decimal separator
  if (point_count == 0 && comma_count > 1) {
    return true;
  }

  return false;
}

bool StringToDecimalConverter::is_comma_as_separator(
  const std::vector<std::string> &strs)
{
  std::vector<std::string> comma_as_separator_strs;
  std::vector<std::string> point_as_separator_strs;

  for (const auto &str : strs) {
    bool is_comma_sep_for_sure = is_comma_as_decimal_separator(str);
    bool is_point_sep_for_sure = is_point_as_decimal_separator(str);

    if (is_comma_sep_for_sure) {
      comma_as_separator_strs.push_back(str);
    }
    if (is_point_sep_for_sure) {
      point_as_separator_strs.push_back(str);
    }
  }

  if (comma_as_separator_strs.empty()) {
    return false;
  }

  if (!point_as_separator_strs.empty()) {
    std::cerr << "WARNING: Cannot determine separator with certainty.\n";
    std::cerr << "Comma as separator areas:\n";
    std::cerr << "Contradictory example areas:\n";
    std::cerr << "  Comma as separator: " << comma_as_separator_strs[0]
              << "\n";
    std::cerr << "  Point as separator: " << point_as_separator_strs[0]
              << "\n";
    return false;
  }

  assert(!comma_as_separator_strs.empty());
  assert(!point_as_separator_strs.empty());

  return true;
}

std::string StringToDecimalConverter::parse_str(
  const std::string &str,
  bool is_point_as_separator)
{
  if (is_str_NA(str)) {
    return "-1.0";
  }

  assert(is_str_correct_format(str));

  std::string processed_str = str;

  if (is_point_as_separator) {
    processed_str = remove_char(processed_str, comma_);

    if (count_char(processed_str, point_) > 1) {
      size_t last_point_pos = processed_str.rfind(point_);
      std::string cleaned_str;

      // Only keep the last point. This should not arise in practice if data is
      // consistent
      for (size_t i = 0; i < processed_str.size(); ++i) {
        if (processed_str[i] != point_ || i == last_point_pos) {
          cleaned_str += processed_str[i];
        }
      }

      processed_str = cleaned_str;
    }

  } else {
    processed_str = remove_char(processed_str, point_);
    processed_str[processed_str.rfind(comma_)] = point_;
    processed_str = remove_char(processed_str, comma_);
  }

  return processed_str;
}
