#include "constants.h"
#include "inset_state.h"

// Functions to project map with the Smyth equal-surface projection (also
// known as Craster rectangular projection):
// https://en.wikipedia.org/wiki/Cylindrical_equal-area_projection
// The purpose is to create a projection that has a 2:1 aspect ratio so that
// the Fourier transforms work optimally when padding is reduced to zero.
Point point_after_smyth_craster_projection(Point p1)
{
  return Point(
    p1.x() * std::sqrt(2.0 * pi) / 180.0,
    std::sin(p1.y() * pi / 180.0) * std::sqrt(0.5 * pi));
}

void InsetState::apply_smyth_craster_projection()
{
  std::cerr << "Applying Smyth-Craster projection." << std::endl;
  transform_points(point_after_smyth_craster_projection);
  return;
}

// Functions for the projection from Smyth-Craster coordinates to longitude
// latitude. We assume that the Smyth-Craster coordinates have been scaled
// to fit in the box [0, lx] * [0, ly].
Point point_before_smyth_craster_projection(
  Point p1,
  const unsigned int lx,
  const unsigned int ly)
{
  return Point(
    (p1.x() * 360.0 / lx) - 180.0,
    180.0 * std::asin((2.0 * p1.y() / ly) - 1) / pi);
}

void InsetState::revert_smyth_craster_projection()
{
  // Specialise point_before_smyth_craster_projection with lx_ and ly_
  std::function<Point(Point)> lambda = [lx = lx_, ly = ly_](Point p1) {
    return point_before_smyth_craster_projection(p1, lx, ly);
  };

  // Apply `lambda` to all points
  transform_points(lambda);
  return;
}
