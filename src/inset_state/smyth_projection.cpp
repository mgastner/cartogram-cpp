#include "constants.hpp"
#include "inset_state.hpp"

// Functions to project map with the Smyth equal-surface projection (also
// known as Craster rectangular projection):
// https://en.wikipedia.org/wiki/Cylindrical_equal-area_projection
// The purpose is to create a projection that has a 2:1 aspect ratio so that
// the Fourier transforms work optimally when padding is reduced to zero.
static Point point_after_smyth_craster_projection(const Point &p1)
{
  return Point(
    p1.x() * std::sqrt(2.0 * pi) / 180.0,
    std::sin(p1.y() * pi / 180.0) * std::sqrt(0.5 * pi));
}

void InsetState::apply_smyth_craster_projection()
{
  std::cerr << "Applying Smyth-Craster projection." << std::endl;
  transform_points(point_after_smyth_craster_projection);
}
