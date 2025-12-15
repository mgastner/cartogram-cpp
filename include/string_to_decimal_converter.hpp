#ifndef STRING_TO_DECIMAL_CONVERTER_H
#define STRING_TO_DECIMAL_CONVERTER_H

#include <string>
#include <vector>

/**
 * @brief A utility class for converting string representations of numbers to decimal format.
 * 
 * This class handles various number formats including:
 * - Regular decimal numbers (e.g. "123.456", "123,456")
 * - Scientific notation (e.g. "1.23e-4", "1.23E4")
 * - Special value "NA"
 * 
 * For scientific notation:
 * - Both 'e' and 'E' are supported as exponent markers
 * - The mantissa can use either '.' or ',' as decimal separator
 * - The exponent must be an integer and can be negative
 * - The mantissa must be a valid decimal number
 */
class StringToDecimalConverter
{
private:
  static constexpr char point_ = '.';
  static constexpr char comma_ = ',';
  static constexpr char minus_ = '-';
  static const std::string NA_;

  static bool is_valid_char(char ch);
  static std::string remove_char(std::string str, char ch);
  static bool has_multiple_commas_and_points(const std::string &str);
  static bool has_separator_at_the_end(const std::string &str);
  static bool has_invalid_comma_point_sequence(const std::string &str);

public:
  static bool is_comma_as_separator(const std::vector<std::string> &strs);
  static bool is_str_NA(const std::string &str);
  static bool is_str_valid_characters(const std::string &str);
  static bool is_str_correct_format(const std::string &str);
  static std::string parse_str(const std::string &str, bool);
};

#endif  // STRING_TO_DECIMAL_CONVERTER_H
