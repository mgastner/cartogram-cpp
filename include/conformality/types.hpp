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

inline constexpr Mat2 identity_mat2() noexcept
{
  return {1.0, 0.0, 0.0, 1.0};
}

inline constexpr double det(const Mat2 &a) noexcept
{
  return a.a11 * a.a22 - a.a12 * a.a21;
}

// Squared Frobenius norm
inline constexpr double frobenius_sq(const Mat2 &a) noexcept
{
  return a.a11 * a.a11 + a.a12 * a.a12 + a.a21 * a.a21 + a.a22 * a.a22;
}

inline constexpr Mat2 operator*(const Mat2 &a, const Mat2 &b) noexcept
{
  return {
    a.a11 * b.a11 + a.a12 * b.a21,
    a.a11 * b.a12 + a.a12 * b.a22,
    a.a21 * b.a11 + a.a22 * b.a21,
    a.a21 * b.a12 + a.a22 * b.a22};
}

}  // namespace conformality

#endif  // CONFORMALITY_TYPES_HPP_
