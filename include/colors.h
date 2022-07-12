#ifndef COLORS_H_
#define COLORS_H_

#include <cstdint>
#include <string>
#include <unordered_map>

struct Color {

  // red, green, blue values out of 255
  int r;
  int g;
  int b;
  Color();
  Color(int, int, int);
  explicit Color(std::string);
  std::string eps() const;

  bool operator==(const Color &rhs) const
  {
    return (r == rhs.r && g == rhs.g && b == rhs.b);
  }
};

#endif
