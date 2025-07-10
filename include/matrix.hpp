#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include "cgal_typedef.hpp"

struct Matrix {
  double p11, p12, p13, p21, p22, p23, p31, p32, p33;  // Matrix elements
  Matrix();  // Constructor for identity matrix

  // Convert triangle to matrix
  Matrix(const Point &, const Point &, const Point &);
  void scale(const double);
  double det() const;
  Matrix adjugate() const;
  Matrix inverse() const;
  Matrix multiplied_with(const Matrix &) const;
  Point transformed_point(const Point &) const;
};

#endif  // MATRIX_HPP_
