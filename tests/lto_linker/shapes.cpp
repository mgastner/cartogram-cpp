#include "shapes.hpp"

namespace
{
struct Square final : Shape {
  int sides() const override
  {
    return 4;
  }
};
}  // namespace

std::unique_ptr<Shape> make_square()
{
  return std::make_unique<Square>();
}
