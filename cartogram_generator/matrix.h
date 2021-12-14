#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include "inset_state.h"

struct Matrix {
  double p11, p12, p13, p21, p22, p23, p31, p32, p33;  // Matrix elements
  Matrix();  // Constructor for identity matrix

  // Convert triangle to matrix
  Matrix(const Point a, const Point b, const Point c);
  void scale(const double);
  double det();
  Matrix adjugate();
  Matrix inverse();
  Matrix multiplied_with(const Matrix);
  Point transformed_point(const Point);
};

#endif
