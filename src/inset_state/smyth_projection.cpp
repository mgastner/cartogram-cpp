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

void InsetState::snap_vertices_across_antimeridian()
{
  auto on_meridian = [](const Point &p, const double lon) {
    return !(p.x() < lon) && !(p.x() > lon);
  };
  std::vector<double> west_latitudes;
  auto collect = [&](const Polygon &ring) {
    for (const Point &p : ring) {
      if (on_meridian(p, -180.0))
        west_latitudes.push_back(p.y());
    }
  };
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      collect(pwh.outer_boundary());
      for (const auto &h : pwh.holes())
        collect(h);
    }
  }
  std::sort(west_latitudes.begin(), west_latitudes.end());
  unsigned int n_snapped = 0;
  transform_points([&](const Point &p) {
    if (!on_meridian(p, 180.0) || west_latitudes.empty())
      return p;
    const auto it =
      std::lower_bound(west_latitudes.begin(), west_latitudes.end(), p.y());
    double nearest =
      (it == west_latitudes.end()) ? west_latitudes.back() : *it;
    if (
      it != west_latitudes.begin() &&
      std::fabs(*std::prev(it) - p.y()) < std::fabs(nearest - p.y()))
      nearest = *std::prev(it);
    const double gap = std::fabs(nearest - p.y());
    if (gap > 0.0 && gap < antimeridian_snap_tolerance) {
      ++n_snapped;
      return Point(p.x(), nearest);
    }
    return p;
  });
  if (n_snapped > 0) {
    std::cerr << "Moved " << n_snapped
              << " vertex occurrence(s) at 180 degrees east onto the "
              << "matching vertex at 180 degrees west" << std::endl;
  }
}

void InsetState::apply_smyth_craster_projection()
{
  std::cerr << "Applying Smyth-Craster projection." << std::endl;
  transform_points(point_after_smyth_craster_projection);
}
