#include "constants.hpp"
#include "inset_state.hpp"

// Functions to project map with the Smyth equal-surface projection (also
// known as Craster rectangular projection):
// https://en.wikipedia.org/wiki/Cylindrical_equal-area_projection
// The purpose is to create a projection that has a 2:1 aspect ratio so that
// the Fourier transforms work optimally when padding is reduced to zero.
Point point_after_smyth_craster_projection(const Point &p1)
{
  return Point(
    p1.x() * std::sqrt(2.0 * pi) / 180.0,
    std::sin(p1.y() * pi / 180.0) * std::sqrt(0.5 * pi));
}

std::string InsetState::geo_div_with_edge_across_antimeridian() const
{
  auto crosses = [](const Polygon &ring) {
    for (std::size_t k = 0; k < ring.size(); ++k) {
      const Point &a = ring[k];
      const Point &b = ring[(k + 1 == ring.size()) ? 0 : k + 1];
      if (std::fabs(a.x() - b.x()) > 180.0)
        return true;
    }
    return false;
  };
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      if (crosses(pwh.outer_boundary()))
        return gd.id();
      for (const auto &h : pwh.holes()) {
        if (crosses(h))
          return gd.id();
      }
    }
  }
  return "";
}

void InsetState::apply_smyth_craster_projection()
{
  std::cerr << "Applying Smyth-Craster projection." << std::endl;
  transform_points(point_after_smyth_craster_projection);
}
