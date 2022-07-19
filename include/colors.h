#ifndef COLORS_H_
#define COLORS_H_

#include <cstdint>
#include <string>
#include <unordered_map>

struct Color {

  // red, green, blue values between 0 and 1
  double r;
  double g;
  double b;

  // Defaults to white
  Color()
  {
    r = 1.0;
    g = 1.0;
    b = 1.0;
  };
  Color(double red, double green, double blue);
  explicit Color(std::string color_as_string);

  bool operator==(const Color &rhs) const
  {
    return (r == rhs.r && g == rhs.g && b == rhs.b);
  }

  // Getter
  double operator()(const char c) const
  {
    switch (c) {
    case 'r':
      return r;
    case 'g':
      return g;
    case 'b':
      return b;

    // Return -1.0 (incorrect color)
    default:
      return -1.0;
    }
  }

  // Setter
  double &operator()(const char c)
  {
    switch (c) {
    case 'r':
      return r;
    case 'g':
      return g;
    case 'b':
      return b;

    // Return r by default
    default:
      return r;
    }
  }
};

#endif
