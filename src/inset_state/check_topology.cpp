#include "inset_state.hpp"

// Returns error if there are holes not inside their respective polygons
void InsetState::holes_inside_polygons() const
{
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto &ext_ring = pwh.outer_boundary();
      for (auto &h : pwh.holes()) {
        for (auto &p : h) {

          // TODO: In the future, a better method would be to only check
          // whether one point in each hole is on the bounded side of the
          // exterior ring. Next, check whether the hole intersects the
          // polygon at any point. For this, the function
          // "do_intersect(Polygon, Polygon)" may help.
          if (ext_ring.bounded_side(p) == CGAL::ON_UNBOUNDED_SIDE) {
            CGAL::set_pretty_mode(std::cerr);
            std::cerr << "Hole detected outside polygon!" << std::endl;
            std::cerr << "Hole: " << h << std::endl;
            std::cerr << "Polygon: " << ext_ring << std::endl;
            std::cerr << "GeoDiv: " << gd.id() << std::endl;
            _Exit(20);
          }
        }
      }
    }
  }
}

void InsetState::rings_are_simple() const
{
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto &ext_ring = pwh.outer_boundary();
      if (!ext_ring.is_simple()) {
        std::cerr << "External ring not a simple polygon!" << std::endl;
        std::cerr << "Coordinates: " << ext_ring << std::endl;
        std::cerr << "GeoDiv: " << gd.id() << std::endl;
        _Exit(43567);
      }
      for (const auto &h : pwh.holes()) {
        if (!h.is_simple()) {
          std::cerr << "Hole is not a simple polygon!" << std::endl;
          std::cerr << "Coordinates: " << h << std::endl;
          std::cerr << "GeoDiv: " << gd.id() << std::endl;
          _Exit(43568);
        }
      }
    }
  }
}

void InsetState::check_topology() const
{
  holes_inside_polygons();
  rings_are_simple();
}
