#ifndef COLORS_H_
#define COLORS_H_

#include <cstdint>
#include <unordered_map>
#include <string>

struct Color {

  // red, green, blue values out of 255
  int r;
  int g;
  int b;
  Color();
  Color(int, int, int);
  Color(std::string);
  std::string eps();

  bool operator == (const Color& rhs)
  {
      return (r == rhs.r && g == rhs.g && b == rhs.b);
  }
};

static std::unordered_map<std::string, std::string> html_colors =
{

  // From https://en.wikipedia.org/wiki/Web_colors
  // Shades of pink
  {"mediumvioletred", "#C71585"},
  {"deeppink", "#FF1493"},
  {"palevioletred", "#DB7093"},
  {"hotpink", "#FF69B4"},
  {"lightpink", "#FFB6C1"},
  {"pink", "#FFC0CB"},

  // Shades of red
  {"darkred", "#8B0000"},
  {"red", "#FF0000"},
  {"firebrick", "#B22222"},
  {"crimson", "#DC143C"},
  {"indianred", "#CD5C5C"},
  {"lightcoral", "#F08080"},
  {"salmon", "#FA8072"},
  {"darksalmon", "#E9967A"},
  {"lightsalmon", "#FFA07A"},

  // Shades of orange
  {"orangered", "#FF4500"},
  {"tomato", "#FF6347"},
  {"darkorange", "#FF8C00"},
  {"coral", "#FF7F50"},
  {"orange", "#FFA500"},

  // Shades of yellow
  {"darkkhaki", "#BDB76B"},
  {"gold", "#FFD700"},
  {"khaki", "#F0E68C"},
  {"peachpuff", "#FFDAB9"},
  {"yellow", "#FFFF00"},
  {"palegoldenrod", "#EEE8AA"},
  {"moccasin", "#FFE4B5"},
  {"papayawhip", "#FFEFD5"},
  {"lightgoldenrodyellow", "#FAFAD2"},
  {"lemonchiffon", "#FFFACD"},
  {"lightyellow", "#FFFFE0"},

  // Shades of brown
  {"maroon", "#800000"},
  {"brown", "#A52A2A"},
  {"saddlebrown", "#8B4513"},
  {"sienna", "#A0522D"},
  {"chocolate", "#D2691E"},
  {"darkgoldenrod", "#B8860B"},
  {"peru", "#CD853F"},
  {"rosybrown", "#BC8F8F"},
  {"goldenrod", "#DAA520"},
  {"sandybrown", "#F4A460"},
  {"tan", "#D2B48C"},
  {"burlywood", "#DEB887"},
  {"wheat", "#F5DEB3"},
  {"navajowhite", "#FFDEAD"},
  {"bisque", "#FFE4C4"},
  {"blanchedalmond", "#FFEBCD"},
  {"cornsilk", "#FFF8DC"},

  // Shades of green
  {"darkgreen", "#006400"},
  {"green", "#008000"},
  {"darkolivegreen", "#556B2F"},
  {"forestgreen", "#228B22"},
  {"seagreen", "#2E8B57"},
  {"olive", "#808000"},
  {"olivedrab", "#6B8E23"},
  {"mediumseagreen", "#3CB371"},
  {"limegreen", "#32CD32"},
  {"lime", "#00FF00"},
  {"springgreen", "#00FF7F"},
  {"mediumspringgreen", "#00FA9A"},
  {"darkseagreen", "#8FBC8F"},
  {"mediumaquamarine", "#66CDAA"},
  {"yellowgreen", "#9ACD32"},
  {"lawngreen", "#7CFC00"},
  {"chartreuse", "#7FFF00"},
  {"lightgreen", "#90EE90"},
  {"greenyellow", "#ADFF2F"},
  {"palegreen", "#98FB98"},

  // Shades of cyan
  {"teal", "#008080"},
  {"darkcyan", "#008B8B"},
  {"lightseagreen", "#20B2AA"},
  {"cadetblue", "#5F9EA0"},
  {"darkturquoise", "#00CED1"},
  {"mediumturquoise", "#48D1CC"},
  {"turquoise", "#40E0D0"},
  {"aqua", "#00FFFF"},
  {"cyan", "#00FFFF"},
  {"aquamarine", "#7FFFD4"},
  {"paleturquoise", "#AFEEEE"},
  {"lightcyan", "#E0FFFF"},

  // Shades of blue
  {"navy", "#000080"},
  {"darkblue", "#00008B"},
  {"mediumblue", "#0000CD"},
  {"blue", "#0000FF"},
  {"midnightblue", "#191970"},
  {"royalblue", "#4169E1"},
  {"steelblue", "#4682B4"},
  {"dodgerblue", "#1E90FF"},
  {"deepskyblue", "#00BFFF"},
  {"cornflowerblue", "#6495ED"},
  {"skyblue", "#87CEEB"},
  {"lightskyblue", "#87CEFA"},
  {"lightsteelblue", "#B0C4DE"},
  {"lightblue", "#ADD8E6"},
  {"powderblue", "#B0E0E6"},

  // Shades of purple, violet, and magenta
  {"indigo", "#4B0082"},
  {"purple", "#800080"},
  {"darkmagenta", "#8B008B"},
  {"darkviolet", "#9400D3"},
  {"darkslateblue", "#483D8B"},
  {"blueviolet", "#8A2BE2"},
  {"darkorchid", "#9932CC"},
  {"fuchsia", "#FF00FF"},
  {"magenta", "#FF00FF"},
  {"slateblue", "#6A5ACD"},
  {"mediumslateblue", "#7B68EE"},
  {"mediumorchid", "#BA55D3"},
  {"mediumpurple", "#9370DB"},
  {"orchid", "#DA70D6"},
  {"violet", "#EE82EE"},
  {"plum", "#DDA0DD"},
  {"thistle", "#D8BFD8"},
  {"lavender", "#E6E6FA"},

  // Shades of white
  {"mistyrose", "#FFE4E1"},
  {"antiquewhite", "#FAEBD7"},
  {"linen", "#FAF0E6"},
  {"beige", "#F5F5DC"},
  {"whitesmoke", "#F5F5F5"},
  {"lavenderblush", "#FFF0F5"},
  {"oldlace", "#FDF5E6"},
  {"aliceblue", "#F0F8FF"},
  {"seashell", "#FFF5EE"},
  {"ghostwhite", "#F8F8FF"},
  {"honeydew", "#F0FFF0"},
  {"floralwhite", "#FFFAF0"},
  {"azure", "#F0FFFF"},
  {"mintcream", "#F5FFFA"},
  {"snow", "#FFFAFA"},
  {"ivory", "#FFFFF0"},
  {"white", "#FFFFFF"},

  // Shades of gray and black
  {"black", "#000000"},
  {"darkslategray", "#2F4F4F"},
  {"dimgray", "#696969"},
  {"slategray", "#708090"},
  {"gray", "#808080"},
  {"lightslategray", "#778899"},
  {"darkgray", "#A9A9A9"},
  {"silver", "#C0C0C0"},
  {"lightgray", "#D3D3D3"},
  {"gainsboro", "#DCDCDC"}
};

#endif
