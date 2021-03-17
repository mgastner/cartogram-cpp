#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include "map_state.h"

struct Matrix {

  // all positions of matrix
  double p11, p12, p13, p21, p22, p23, p31, p32, p33;
  Matrix(); // constructor for identity matrix
  Matrix (const Matrix &); // copy constructor
  Matrix(Triangle); // matrix from a triangle
  Matrix(XYPoint, XYPoint, XYPoint); // matrix from triangle of XYPoints
  void scale(double); //
  void print();
  double det();
  Matrix adjugate();
  Matrix inverse();
  Matrix multiply(Matrix);
  Point transform_point(Point);
  XYPoint transform_XYPoint(XYPoint);

};

#endif
