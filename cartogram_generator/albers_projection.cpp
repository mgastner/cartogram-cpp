#include <iostream>
#include <utility>

#include "cgal_typedef.h"
#include "map_state.h"

void albers_projection(InsetState *inset_state) {
  std::vector<GeoDiv> inset_state_converted = inset_state->geo_divs();

  for (auto gd : inset_state_converted) {
    for (auto pwh : gd.polygons_with_holes()) {
      Polygon pwh_outer = pwh.outer_boundary();

      for (unsigned int i = 0; i < pwh_outer.size(); i++) {
      }

      for (auto hole_it = pwh.holes_begin(); hole_it != pwh.holes_end(); hole_it++) {
        auto hole = *hole_it;
        for (unsigned int i = 0; i < hole.size(); i++) {
        }
      }
    }
  }

  inset_state->set_geo_divs(inset_state_converted);
}