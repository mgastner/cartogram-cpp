#include "cartogram_info.h"
#include "inset_state.h"
#include "round_point.h"
#include <CGAL/Boolean_set_operations_2.h>

// Returns error if there are holes not inside their respective polygons
void holes_inside_polygons(InsetState *inset_state)
{
  for (const auto &gd : inset_state->geo_divs()) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto ext_ring = pwh.outer_boundary();
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        for (unsigned int i = 0; i < h->size(); ++i) {

          // TODO: In the future, a better method would be to only check
          // whether one point in each hole is on the bounded side of the
          // exterior ring. Next, check whether the hole intersects the
          // polygon at any point. For this, the function
          // "do_intersect(Polygon, Polygon)" may help.
          if (ext_ring.bounded_side((*h)[i]) == CGAL::ON_UNBOUNDED_SIDE) {
            CGAL::set_pretty_mode(std::cerr);
            std::cerr << "Hole detected outside polygon!" << std::endl;
            std::cerr << "Hole: " << (*h) << std::endl;
            std::cerr << "Polygon: " << ext_ring << std::endl;
            std::cerr << "GeoDiv: " << gd.id() << std::endl;
            _Exit(20);
          }
        }
      }
    }
  }
  return;
}

bool duplicates(std::vector<Point> v) {
  CGAL::set_pretty_mode(std::cerr);
  for (size_t i = 0; i < v.size() - 1; ++i) {
    if (points_almost_equal(v[i], v[i + 1])) {
      std::cerr << "i = " << i << std::endl;
      std::cerr << "Point: " << i << ", v[i]: " << v[i] << std::endl;
      std::cerr << "Point: "
                << i + 1
                << ", v[i + 1]: "
                << v[i + 1]
                << std::endl;
      return true;
    }
  }
  return false;
}

void rings_are_simple(InsetState *inset_state)
{
  for (const auto &gd : inset_state->geo_divs()) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto ext_ring = pwh.outer_boundary();
      if (!ext_ring.is_simple()) {
        std::cerr << "External ring not a simple polygon!" << std::endl;
        std::cerr << "Coordinates: " << ext_ring << std::endl;
        std::cerr << "GeoDiv: " << gd.id() << std::endl;
        _Exit(43567);
      }
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        if (!h->is_simple()) {
          std::cerr << "Hole is not a simple polygon!" << std::endl;
          std::cerr << "Coordinates: " << (*h) << std::endl;
          std::cerr << "GeoDiv: " << gd.id() << std::endl;
          _Exit(43568);
        }
      }
    }
  }
  return;
}
