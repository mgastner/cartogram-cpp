#include "interpolate_bilinearly.hpp"
#include <iostream>

// TODO: REPLACE WITH LINEAR INTERPOLATION BASED ON TRIANGULATION

// Function to bilinearly interpolate a numerical array
// grid[0..lx-1][0..*ly-1] whose entries are numbers for the positions:
// x = (0.5, 1.5, ..., lx-0.5), y = (0.5, 1.5, ..., ly-0.5).
// The final argument "zero" can take two possible values: 'x' or 'y'. If
// zero == x, the interpolated function is forced to return 0 if x=0 or x=lx.
// This option is suitable fo interpolating from gridvx because there can be
// no flow through the boundary. If zero == y, the interpolation returns 0 if
// y=0 or y=ly, suitable for gridvy. The unconstrained boundary will be
// determined by continuing the function value at 0.5 (or lx-0.5 or ly-0.5)
// all the way to the edge (i.e. the slope is 0 consistent with a cosine
// transform).
double interpolate_bilinearly(
  const double x,
  const double y,
  const boost::multi_array<double, 2> &grid,
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
    fx0y0 = grid[static_cast<int>(x0)][static_cast<int>(y0)];
  }

  // Function value at (x0, y1).
  double fx0y1;
  if (
    (x < 0.5 && y >= ly - 0.5) || (x < 0.5 && zero == 'x') ||
    (y >= ly - 0.5 && zero == 'y')) {
    fx0y1 = 0.0;
  } else if (x >= 0.5 && y >= ly - 0.5 && zero == 'x') {
    fx0y1 = grid[static_cast<int>(x0)][ly - 1];
  } else {
    fx0y1 = grid[static_cast<int>(x0)][static_cast<int>(y1)];
  }

  // Function value at (x1, y0).
  double fx1y0;
  if (
    (x >= lx - 0.5 && y < 0.5) || (x >= lx - 0.5 && zero == 'x') ||
    (y < 0.5 && zero == 'y')) {
    fx1y0 = 0.0;
  } else if (x >= lx - 0.5 && y >= 0.5 && zero == 'y') {
    fx1y0 = grid[lx - 1][static_cast<int>(y0)];
  } else {
    fx1y0 = grid[static_cast<int>(x1)][static_cast<int>(y0)];
  }

  // Function value at (x1, y1).
  double fx1y1;
  if (
    (x >= lx - 0.5 && y >= ly - 0.5) || (x >= lx - 0.5 && zero == 'x') ||
    (y >= ly - 0.5 && zero == 'y')) {
    fx1y1 = 0.0;
  } else if (x >= lx - 0.5 && y < ly - 0.5 && zero == 'y') {
    fx1y1 = grid[lx - 1][static_cast<int>(y1)];
  } else if (x < lx - 0.5 && y >= ly - 0.5 && zero == 'x') {
    fx1y1 = grid[static_cast<int>(x1)][ly - 1];
  } else {
    fx1y1 = grid[static_cast<int>(x1)][static_cast<int>(y1)];
  }
  return (1.0 - delta_x) * (1.0 - delta_y) * fx0y0 +
         (1.0 - delta_x) * delta_y * fx0y1 +
         delta_x * (1.0 - delta_y) * fx1y0 + delta_x * delta_y * fx1y1;
}

// Overloaded function for flatten_density_with_node_vertices().
double interpolate_bilinearly(
  const double x,
  const double y,
  std::function<double(unsigned int, unsigned int, char)> &cal_velocity,
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
