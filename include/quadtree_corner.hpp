#pragma once

#include "cgal_typedef.hpp"

struct QuadtreeCorner {
public:
  QuadtreeCorner(uint32_t x, uint32_t y) : x_{x}, y_{y} {}

  [[nodiscard]] uint32_t x() const noexcept
  {
    return x_;
  }
  [[nodiscard]] uint32_t y() const noexcept
  {
    return y_;
  }

  [[nodiscard]] operator Point() const noexcept
  {
    return Point(static_cast<double>(x_), static_cast<double>(y_));
  }

  auto operator<=>(const QuadtreeCorner &) const noexcept = default;
  bool operator==(const QuadtreeCorner &) const noexcept = default;

private:
  uint32_t x_;
  uint32_t y_;
};

inline std::ostream &operator<<(std::ostream &os, QuadtreeCorner const &c)
{
  return os << "QuadtreeCorner(" << c.x() << ',' << c.y() << ')';
}
