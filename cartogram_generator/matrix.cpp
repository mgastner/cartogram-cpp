#include "matrix.h"
#include "map_state.h"


Matrix::Matrix() {

  // identity matrix
  p11 = 1;
  p12 = 0;
  p13 = 0;
  p21 = 0;
  p22 = 1;
  p23 = 0;
  p31 = 0;
  p32 = 0;
  p33 = 1;
}

// Matrix from a triangle
// Matrix::Matrix (Triangle triangle) {
//
//   // first vertex
//   p11 = triangle[0][0];
//   p21 = triangle[0][1];
//
//   // second vertex
//   p12 = triangle[1][0];
//   p22 = triangle[1][1];
//
//   // third vertex
//   p13 = triangle[2][0];
//   p23 = triangle[2][1];
//
//   // to make it a 3x3 matrix
//   p31 = 1;
//   p32 = 1;
//   p33 = 1;
// }

// Matrix from three XYPoints
Matrix::Matrix (XYPoint a, XYPoint b, XYPoint c) {

  // first vertex
  p11 = a.x;
  p21 = a.y;

  // second vertex
  p12 = b.x;
  p22 = b.y;

  // third vertex
  p13 = c.x;
  p23 = c.y;

  // to make it a 3x3 matrix
  p31 = 1;
  p32 = 1;
  p33 = 1;
}

// Matrix from vector of three XYPoints

void Matrix::scale(double multiplier) {

  p11 *= multiplier;
  p12 *= multiplier;
  p13 *= multiplier;
  p21 *= multiplier;
  p22 *= multiplier;
  p23 *= multiplier;
  p31 *= multiplier;
  p32 *= multiplier;
  p33 *= multiplier;
}

// For debugging
void Matrix::print() {
  std::cout << p11 << " " << p12 << " " << p13 << "\n\n";
  std::cout << p21 << " " << p22 << " " << p23 << "\n\n";
  std::cout << p31 << " " << p32 << " " << p33 << "\n\n";
}

// calculate determinant
double Matrix::det() {
  return p11 * ((p22 * p33) - (p23 * p32)) -
         p12 * ((p21 * p33) - (p23 * p31)) +
         p13 * ((p21 * p32) - (p22 * p31));
}

Matrix Matrix::adjugate() {
  Matrix adj;

  adj.p11 = ((p22 * p33) - (p23 * p32));
  adj.p21 = - ((p21 * p33) - (p23 * p31));
  adj.p31 = ((p21 * p32) - (p22 * p31));
  adj.p12 = - ((p12 * p33) - (p13 * p32));
  adj.p22 = ((p11 * p33) - (p13 * p31));
  adj.p32 = - ((p11 * p32) - (p12 * p31));
  adj.p13 = ((p12 * p23) - (p13 * p22));
  adj.p23 = - ((p11 * p23) - (p13 * p21));
  adj.p33 = ((p11 * p22) - (p12 * p21));

  return adj;
}

Matrix Matrix::inverse() {

  // calculating adjugate
  Matrix inv = adjugate();

  // dividing by determinant
  inv.scale(1 / det());

  // returning resultant matrix
  return inv;
}

Matrix Matrix::multiply(Matrix m1) {

  Matrix result;

  result.p11 = (p11 * m1.p11) + (p12 * m1.p21) + (p13 * m1.p31);
  result.p12 = (p11 * m1.p12) + (p12 * m1.p22) + (p13 * m1.p32);
  result.p13 = (p11 * m1.p13) + (p12 * m1.p23) + (p13 * m1.p33);
  result.p21 = (p21 * m1.p11) + (p22 * m1.p21) + (p23 * m1.p31);
  result.p22 = (p21 * m1.p12) + (p22 * m1.p22) + (p23 * m1.p32);
  result.p23 = (p21 * m1.p13) + (p22 * m1.p23) + (p23 * m1.p33);
  result.p31 = (p31 * m1.p11) + (p32 * m1.p21) + (p33 * m1.p31);
  result.p32 = (p31 * m1.p12) + (p32 * m1.p22) + (p33 * m1.p32);
  result.p33 = (p31 * m1.p13) + (p32 * m1.p22) + (p33 * m1.p33);

  return result;

}

// transforms point based on a transformation matrix
Point Matrix::transform_point(Point point) {

  double x = p11 * (point[0]) + p12 * (point[1]) + p13;
  double y = p21 * (point[0]) + p22 * (point[1]) + p23;

  Point transformed(x, y);
  return transformed;
}

// transforms XYPoint based on a transformation matrix
XYPoint Matrix::transform_XYPoint(XYPoint point) {

  XYPoint transformed;
  transformed.x = p11 * (point.x) + p12 * (point.y) + p13;
  transformed.y = p21 * (point.x) + p22 * (point.y) + p23;

  return transformed;
}
