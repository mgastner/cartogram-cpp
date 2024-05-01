#include "string_to_decimal_converter.hpp"
#include <cassert>
#include <iostream>

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

  iv) Point appears first (123.456.789,123):
    We assume points are used as big-number separators and remove them
    (123.456.789,123 -> 123456789.123)

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

// Implement the parse_str function with embedded error handling
double StringToDecimalConverter::parse_str(const std::string &area_as_str)
{
    // Check if string is "NA" or empty
    if (area_as_str.empty() || area_as_str == "NA") {
      return -1.0;
    }

    // Check for valid characters in the string
    if (area_as_str.find_first_not_of("0123456789.,-") != std::string::npos) {
      std::cerr << "ERROR: Invalid area string: " << area_as_str << std::endl;
      std::cerr << "Area string must only contain 0-9, '.', '-' and ',' or 'NA'." << std::endl;
      std::exit(18);
    }

    // Validate format regarding commas and points
    size_t n_comma = std::count(area_as_str.begin(), area_as_str.end(), ',');
    size_t n_point = std::count(area_as_str.begin(), area_as_str.end(), '.');

    // Neither commas nor points
    if (n_comma == 0 && n_point == 0) {
      return std::stod(area_as_str);
    }

    // More than one of each indicates error in formatting
    if (n_comma > 1 && n_point > 1) {
      std::cerr << "ERROR: Invalid area string format: " << area_as_str << std::endl;
      std::exit(19);
    }


    std::string processed_str = area_as_str;

    // Contains both commas and points
    if (n_comma > 0 && n_point > 0) {
      if (area_as_str.find('.') < area_as_str.find(',')) {

        // Points are thousands separators, remove them
        processed_str.erase(std::remove(processed_str.begin(), processed_str.end(), '.'), processed_str.end());
        size_t comma_pos = processed_str.find(',');
        processed_str[comma_pos] = '.';

      } else {
        // Commas are thousands separators, remove them
        processed_str.erase(std::remove(processed_str.begin(), processed_str.end(), ','), processed_str.end());
      }
    } else if (n_comma > 0) {
      // If only one comma, likely a decimal
      if (n_comma == 1) {
        size_t comma_pos = processed_str.find(',');
        processed_str[comma_pos] = '.';
      } else {
        processed_str.erase(std::remove(processed_str.begin(), processed_str.end(), ','), processed_str.end());
      }
    } else if (n_point > 0) {
      // If only one decimal, likely a decimal
      processed_str.erase(std::remove(processed_str.begin(), processed_str.end(), '.'), processed_str.end());
    }
    // else, already handled earlier

    double parsed_value = std::stod(processed_str);
    if (parsed_value < 0.0) {
        std::cerr << "ERROR: Negative area in CSV" << std::endl;
        std::exit(101);
    }

    return parsed_value;
}