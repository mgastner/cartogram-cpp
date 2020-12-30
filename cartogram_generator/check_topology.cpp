#include "map_state.h"
#include <CGAL/Boolean_set_operations_2.h>

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
            std::cerr << "Hole: " << hole << std::endl;
            std::cerr << "Polygon: " << ext_ring << std::endl;
            std::cerr << "GeoDiv: " << gd.id() << std::endl;
            _Exit(20);
          }
        }

        // int intersections = 0;
        // double line_y = hole[0][1];
        // double prev_point[2];
        // prev_point[0] = ext_ring[ext_ring.size()-1][0];
        // prev_point[1] = ext_ring[ext_ring.size()-1][1];
        // for (unsigned int i = 0; i < ext_ring.size(); ++i) {
        //   double curr_point[2];
        //   curr_point[0] = ext_ring[i][0];
        //   curr_point[1] = ext_ring[i][1];
        //
        //   // checking if intersection present
        //   if ((curr_point[1] <= line_y && prev_point[1] >= line_y) ||
        //        (curr_point[1] >= line_y && prev_point[1] <= line_y)) {
        //          intersections++; // found an intersection
        //    }
        //    prev_point[0] = curr_point[0];
        //    prev_point[1] = curr_point[1];
        // }

        // There must be an odd number of intersections and the hole must not
        // intersect with the exterior ring
        // if (intersections % 2 == 1 || CGAL::do_intersect(hole, ext_ring)) {
        //   CGAL::set_pretty_mode(std::cout);
        //   std::cerr << "Hole detected outside polygon!" << std::endl;
        //   std::cerr << "Hole: " << hole << std::endl;
        //   std::cerr << "Polygon: " << ext_ring << std::endl;
        //   std::cerr << "GeoDiv: " << gd.id() << std::endl;
        //   _Exit(20);
        // }

      }
    }
  }
  return;
}
