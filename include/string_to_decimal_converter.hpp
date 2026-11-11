#ifndef STRING_TO_DECIMAL_CONVERTER_H
#define STRING_TO_DECIMAL_CONVERTER_H

#include <string>
#include <vector>

// Functions that check the numbers in a CSV file and rewrite them with a
// point as the decimal separator and without group separators. In the file,
// either a point or a comma may be the decimal separator.
namespace StringToDecimalConverter
{
bool is_comma_as_separator(const std::vector<std::string> &strs);
bool is_str_NA(const std::string &str);
bool is_str_valid_characters(const std::string &str);
bool is_str_correct_format(const std::string &str);
std::string parse_str(const std::string &str, bool);
}  // namespace StringToDecimalConverter

#endif  // STRING_TO_DECIMAL_CONVERTER_H
