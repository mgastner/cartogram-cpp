#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include "inset_state.h"

struct Matrix {
  double p11, p12, p13, p21, p22, p23, p31, p32, p33;  // Matrix elements
  Matrix();  // Constructor for identity matrix
  Matrix(XYPoint a, XYPoint b, XYPoint c);  // Convert triangle to matrix
  void scale(double);
  void print();
  double det();
  Matrix adjugate();
  Matrix inverse();
  Matrix multiplied_with(Matrix);
  XYPoint transformed_XYPoint(XYPoint);
};

#endif
