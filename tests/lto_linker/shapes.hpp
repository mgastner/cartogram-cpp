#ifndef SHAPES_HPP_
#define SHAPES_HPP_

#include <memory>

struct Shape {
  virtual ~Shape() = default;
  virtual int sides() const = 0;
};

std::unique_ptr<Shape> make_square();

#endif  // SHAPES_HPP_
