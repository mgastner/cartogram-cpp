#ifndef TARGET_AREA_PARSER_H
#define TARGET_AREA_PARSER_H

#include <string>

class TargetAreaParser
{
private:
  static constexpr char point_ = '.';
  static constexpr char comma_ = ',';
  static constexpr char minus_ = '-';
  static const std::string NA_;

  static bool is_valid_char(char ch);
  static std::string remove_char(std::string str, char ch);
  static int count_char(const std::string &str, char ch);
  static bool has_multiple_commas_and_points(const std::string &area_str);
  static bool has_separator_at_the_end(const std::string &str);
  static bool has_invalid_comma_point_sequence(const std::string &str);

public:
  static bool is_area_str_NA(const std::string &area_str);
  static bool is_area_str_valid_characters(const std::string &area_str);
  static bool is_area_str_correct_format(const std::string &area_str);
  static double parse_area_str(const std::string &area_str);
};

#endif  // TARGET_AREA_PARSER_H
