#ifndef CONFORMALITY_TYPES_HPP_
#define CONFORMALITY_TYPES_HPP_

namespace conformality
{

// Point or vector in the plane
struct Vec2 {
  double x{0.0};
  double y{0.0};
};

// 2x2 matrix. Rows index output components and columns index input
// coordinates, so that for a map f the Jacobian is a_ij = df_i/dx_j.
struct Mat2 {
  double a11{0.0};
  double a12{0.0};
  double a21{0.0};
  double a22{0.0};
};

}  // namespace conformality

#endif  // CONFORMALITY_TYPES_HPP_
