#include "map_state.h"

// Returns true if every hole is inside its respective polygon
bool holes_inside_polygons(MapState *map_state)
{
  for (auto gd : map_state->geo_divs()) {
    for (auto pwh : gd.polygons_with_holes()) {
      Polygon ext_ring = pwh.outer_boundary();
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
        Polygon hole = *hci;
        for (unsigned int i = 0; i < hole.size(); ++i) {
          if (ext_ring.bounded_side(hole[i]) != CGAL::ON_BOUNDED_SIDE) {
            CGAL::set_pretty_mode(std::cout);
            std::cout << "Hole detected outside polygon!" << std::endl;
            std::cout << "Point: " << hole[i] << std::endl;
            std::cout << "Hole: " << hole << std::endl;
            std::cout << "Polygon: " << ext_ring << std::endl;
            std::cout << "GeoDiv: " << gd.id() << std::endl;
            return false;
          }
        }
      }
    }
  }
  return true;
}
