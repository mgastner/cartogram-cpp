#include "constants.h"
#include "matrix.h"
#include "inset_state.h"

Matrix::Matrix()
{
  // Identity matrix
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

// Matrix from three Points
Matrix::Matrix (Point a, Point b, Point c)
{
  // First vertex
  p11 = a.x();
  p21 = a.y();

  // Second vertex
  p12 = b.x();
  p22 = b.y();

  // Third vertex
  p13 = c.x();
  p23 = c.y();

  // Make it a 3x3 matrix
  p31 = 1;
  p32 = 1;
  p33 = 1;
}

void Matrix::scale(double multiplier)
{
  p11 *= multiplier;
  p12 *= multiplier;
  p13 *= multiplier;
  p21 *= multiplier;
  p22 *= multiplier;
  p23 *= multiplier;
  p31 *= multiplier;
  p32 *= multiplier;
  p33 *= multiplier;
  return;
}

// Determinant
double Matrix::det()
{
  return p11 * ((p22 * p33) - (p23 * p32)) -
         p12 * ((p21 * p33) - (p23 * p31)) +
         p13 * ((p21 * p32) - (p22 * p31));
}

Matrix Matrix::adjugate()
{
  Matrix adj;
  adj.p11 = ((p22 * p33) - (p23 * p32));
  adj.p21 = -((p21 * p33) - (p23 * p31));
  adj.p31 = ((p21 * p32) - (p22 * p31));
  adj.p12 = -((p12 * p33) - (p13 * p32));
  adj.p22 = ((p11 * p33) - (p13 * p31));
  adj.p32 = -((p11 * p32) - (p12 * p31));
  adj.p13 = ((p12 * p23) - (p13 * p22));
  adj.p23 = -((p11 * p23) - (p13 * p21));
  adj.p33 = ((p11 * p22) - (p12 * p21));
  return adj;
}

Matrix Matrix::inverse()
{
  // Calculate adjugate
  Matrix inv = adjugate();

  // Divide by determinant
  if (abs(det() < dbl_epsilon)) {
    std::cerr << "ERROR: Matrix inversion for (nearly) singular input"
              << std::endl;
    exit(1);
  }
  inv.scale(1.0 / det());

  // Return resultant matrix
  return inv;
}

Matrix Matrix::multiplied_with(Matrix m1)
{
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

// Transform Point based on a transformation matrix
Point Matrix::transformed_point(Point point)
{
  // XYPoint transformed;
  // transformed.x = p11 * (point.x) + p12 * (point.y) + p13;
  // transformed.y = p21 * (point.x) + p22 * (point.y) + p23;
  // return transformed;
  return Point(
    p11 * (point.x()) + p12 * (point.y()) + p13,
    p21 * (point.x()) + p22 * (point.y()) + p23
  );
}
