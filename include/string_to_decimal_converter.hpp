#ifndef STRING_TO_DECIMAL_CONVERTER_H
#define STRING_TO_DECIMAL_CONVERTER_H

#include <string>

class StringToDecimalConverter
{
private:
  static constexpr char point_ = '.';
  static constexpr char comma_ = ',';
  static constexpr char minus_ = '-';
  static const std::string NA_;

  static bool is_valid_char(char ch);
  static std::string remove_char(std::string str, char ch);
  static int count_char(const std::string &str, char ch);
  static bool has_multiple_commas_and_points(const std::string &str);
  static bool has_separator_at_the_end(const std::string &str);
  static bool has_invalid_comma_point_sequence(const std::string &str);

public:
  static bool is_str_NA(const std::string &str);
  static bool is_str_valid_characters(const std::string &str);
  static bool is_str_correct_format(const std::string &str);
  static double parse_str(const std::string &str);
};

#endif  // STRING_TO_DECIMAL_CONVERTER_H
