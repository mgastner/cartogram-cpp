#ifndef INTERPOL_HPP_
#define INTERPOL_HPP_

#include <algorithm>
#include <cmath>
#include <iostream>

template <typename VelocityCalculator>
double interpolate_bilinearly(
  const double x,
  const double y,
  VelocityCalculator &&cal_velocity,
  const char zero,
  const unsigned int lx,
  const unsigned int ly)
{
  if (x < 0 || x > lx || y < 0 || y > ly) {
    std::cerr << "ERROR: coordinate outside bounding box in " << __func__
              << "(). " << "x=" << x << ", y=" << y << std::endl;
    exit(1);
  }
  if (zero != 'x' && zero != 'y') {
    std::cerr << "ERROR: unknown argument zero in " << __func__ << "()."
              << std::endl;
    exit(1);
  }

  // x0 is the nearest grid point smaller than x.
  // Exception: if x < 0.5, x0 becomes 0.0.
  const double x0 = std::max(0.0, floor(x + 0.5) - 0.5);

  // x1 is the nearest grid point larger than x.
  // Exception: if x > lx-0.5, x1 becomes lx.
  const double x1 = std::min(static_cast<double>(lx), floor(x + 0.5) + 0.5);

  // Similarly for y
  const double y0 = std::max(0.0, floor(y + 0.5) - 0.5);
  const double y1 = std::min(static_cast<double>(ly), floor(y + 0.5) + 0.5);

  // On a scale from 0 to 1, how far is x (or y) away from x0 (or y0)?
  // 1 means x = x1.
  const double delta_x = (x - x0) / (x1 - x0);
  const double delta_y = (y - y0) / (y1 - y0);

  // Function value at (x0, y0).
  double fx0y0;
  if (
    (x < 0.5 && y < 0.5) || (x < 0.5 && zero == 'x') ||
    (y < 0.5 && zero == 'y')) {
    fx0y0 = 0.0;
  } else {
    fx0y0 = cal_velocity(
      static_cast<unsigned int>(x0),
      static_cast<unsigned int>(y0),
      zero);
  }

  // Function value at (x0, y1).
  double fx0y1;
  if (
    (x < 0.5 && y >= ly - 0.5) || (x < 0.5 && zero == 'x') ||
    (y >= ly - 0.5 && zero == 'y')) {
    fx0y1 = 0.0;
  } else if (x >= 0.5 && y >= ly - 0.5 && zero == 'x') {
    fx0y1 = cal_velocity(static_cast<unsigned int>(x0), ly - 1, zero);
  } else {
    fx0y1 = cal_velocity(
      static_cast<unsigned int>(x0),
      static_cast<unsigned int>(y1),
      zero);
  }

  // Function value at (x1, y0).
  double fx1y0;
  if (
    (x >= lx - 0.5 && y < 0.5) || (x >= lx - 0.5 && zero == 'x') ||
    (y < 0.5 && zero == 'y')) {
    fx1y0 = 0.0;
  } else if (x >= lx - 0.5 && y >= 0.5 && zero == 'y') {
    fx1y0 = cal_velocity(
      static_cast<unsigned int>(lx - 1),
      static_cast<unsigned int>(y0),
      zero);
  } else {
    fx1y0 = cal_velocity(
      static_cast<unsigned int>(x1),
      static_cast<unsigned int>(y0),
      zero);
  }

  // Function value at (x1, y1).
  double fx1y1;
  if (
    (x >= lx - 0.5 && y >= ly - 0.5) || (x >= lx - 0.5 && zero == 'x') ||
    (y >= ly - 0.5 && zero == 'y')) {
    fx1y1 = 0.0;
  } else if (x >= lx - 0.5 && y < ly - 0.5 && zero == 'y') {
    fx1y1 = cal_velocity(
      static_cast<unsigned int>(lx - 1),
      static_cast<unsigned int>(y1),
      zero);
  } else if (x < lx - 0.5 && y >= ly - 0.5 && zero == 'x') {
    fx1y1 = cal_velocity(
      static_cast<unsigned int>(x1),
      static_cast<unsigned int>(ly - 1),
      zero);
  } else {
    fx1y1 = cal_velocity(
      static_cast<unsigned int>(x1),
      static_cast<unsigned int>(y1),
      zero);
  }

  return (1.0 - delta_x) * (1.0 - delta_y) * fx0y0 +
         (1.0 - delta_x) * delta_y * fx0y1 +
         delta_x * (1.0 - delta_y) * fx1y0 + delta_x * delta_y * fx1y1;
}

// Bilinear interpolation of a velocity component that is periodic in x with
// period lx and has reflecting boundaries at y = 0 and y = ly. The grid
// values are at the cell centers (i + 0.5, j + 0.5). Every x is allowed: x is
// reduced modulo lx, and between x = lx - 0.5 and x = lx + 0.5 the
// interpolation uses the values in columns lx - 1 and 0. Near y = 0 and
// y = ly, the y-component interpolates toward 0, and the x-component is
// constant, as in interpolate_bilinearly().
template <typename VelocityCalculator>
double interpolate_bilinearly_periodic_x(
  const double x,
  const double y,
  VelocityCalculator &&cal_velocity,
  const char zero,
  const unsigned int lx,
  const unsigned int ly)
{
  if (!std::isfinite(x) || y < 0 || y > ly) {
    std::cerr << "ERROR: coordinate outside bounding box in " << __func__
              << "(). " << "x=" << x << ", y=" << y << std::endl;
    exit(1);
  }
  if (zero != 'x' && zero != 'y') {
    std::cerr << "ERROR: unknown argument zero in " << __func__ << "()."
              << std::endl;
    exit(1);
  }
  const double dlx = lx;

  // Reduce x to [0, lx). Rounding can turn a tiny negative x into lx.
  double x_red = x - dlx * std::floor(x / dlx);
  if (x_red >= dlx) {
    x_red -= dlx;
  }

  // x0 is the nearest cell center at or left of x_red, which is -0.5 (the
  // center of column lx - 1 shifted by -lx) if x_red < 0.5
  const double x0 = std::floor(x_red + 0.5) - 0.5;
  const double delta_x = x_red - x0;
  const unsigned int i0 = (x0 < 0.0) ? lx - 1 : static_cast<unsigned int>(x0);
  const unsigned int i1 = (i0 + 1 == lx) ? 0 : i0 + 1;

  // y0 and y1 are the nearest cell centers below and above y, or the
  // boundary y = 0 or y = ly
  const double y0 = std::max(0.0, std::floor(y + 0.5) - 0.5);
  const double y1 =
    std::min(static_cast<double>(ly), std::floor(y + 0.5) + 0.5);
  const double delta_y = (y - y0) / (y1 - y0);

  // Value in column i at height yy, which is a cell center or a boundary
  auto value = [&](const unsigned int i, const double yy) {
    if (yy < 0.5 || yy > ly - 0.5) {
      if (zero == 'y') {
        return 0.0;
      }
      return static_cast<double>(
        cal_velocity(i, (yy < 0.5) ? 0u : ly - 1, zero));
    }
    return static_cast<double>(
      cal_velocity(i, static_cast<unsigned int>(yy), zero));
  };
  return (1.0 - delta_x) * (1.0 - delta_y) * value(i0, y0) +
         (1.0 - delta_x) * delta_y * value(i0, y1) +
         delta_x * (1.0 - delta_y) * value(i1, y0) +
         delta_x * delta_y * value(i1, y1);
}

#endif  // INTERPOL_HPP_
