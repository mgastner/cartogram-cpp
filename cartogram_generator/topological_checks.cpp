#include "map_state.h"

// Returns error if there are holes not inside their respective polygons
void holes_inside_polygons(MapState *map_state)
{
  for (auto gd : map_state->geo_divs()) {
    for (auto pwh : gd.polygons_with_holes()) {
      Polygon ext_ring = pwh.outer_boundary();
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
        Polygon hole = *hci;
        for (unsigned int i = 0; i < hole.size(); ++i) {
          if (ext_ring.bounded_side(hole[i]) != CGAL::ON_BOUNDED_SIDE) {
            CGAL::set_pretty_mode(std::cout);
            std::cerr << "Hole detected outside polygon!" << std::endl;
            std::cerr << "Point: " << hole[i] << std::endl;
            std::cerr << "Hole: " << hole << std::endl;
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
