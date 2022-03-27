#include "colors.h"
#include <iostream>
#include <sstream>

Color::Color(int red, int green, int blue) :
  r(red),
  g(green),
  b(blue)
{
  return;
}

Color::Color(std::string color_as_string)
{

  // Make sure colors are case-insensitive
  std::transform(color_as_string.begin(), color_as_string.end(),
                 color_as_string.begin(), ::tolower);
  if (html_colors.find(color_as_string) != html_colors.end()) {

    // HTML color in hexadecimal format
    color_as_string = html_colors[color_as_string];
  }
  if (color_as_string[0] == '#') {

    // Hex code for color
    color_as_string.erase(0, 1);

    // From https://gist.github.com/bert/998020
    int hex_int = stoi(color_as_string, nullptr, 16);
    r = ((hex_int >> 16) & 0xFF); // Extract the RR byte
    g = ((hex_int >> 8) & 0xFF); // Extract the GG byte
    b = ((hex_int) & 0xFF); // Extract the BB byte
  } else if ("rgb" == color_as_string.substr(0, 3)) {

    // RGB value
    std::stringstream css(color_as_string);
    css.ignore(32, '(');
    css >> r;
    if (css.peek() == ',') {
      css.ignore(1, ',');
    }
    css >> g;
    if (css.peek() == ',') {
      css.ignore(1, ',');
    }
    css >> b;
  } else if (isdigit(color_as_string[0])) {

    // We assume the color is given as RGB value
    std::stringstream css(color_as_string);
    css >> r;
    if (css.peek() == ',') {
      css.ignore(1, ',');
    }
    css >> g;
    if (css.peek() == ',') {
      css.ignore(1, ',');
    }
    css >> b;
  } else {

    // Wrong format. We fill the GeoDIV white.
    r = 255;
    g = 255;
    b = 255;
    std::cerr << "WARNING: Color in wrong format. "
              << "Please refer to README.md!"
              << std::endl;
    std::cerr << "Color: " << color_as_string << std::endl;
  }
  return;
}

std::string Color::eps()
{
  std::string temp;
  temp.append(std::to_string(static_cast<double>(r) / 255.0));
  temp.append(" ");
  temp.append(std::to_string(static_cast<double>(g) / 255.0));
  temp.append(" ");
  temp.append(std::to_string(static_cast<double>(b) / 255.0));
  temp.append(" ");
  return temp;
}
