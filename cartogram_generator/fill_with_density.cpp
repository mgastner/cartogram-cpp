#include "map_state.h"
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Bbox_2.h>

void fill_with_density(MapState *map_state)
{
  FTReal2d &rho_init = *map_state->ref_to_rho_init();
  for (unsigned int i = 0; i < map_state->lx(); i++) {
    for (unsigned int j = 0; j < map_state->ly(); j++) {
      rho_init(i, j) = 0;
    }
  }

  std::cout << "Breakpoint 1" << std::endl;

  for (auto gd : map_state->geo_divs()) {

    std::cout << "Working on gd with ID " << gd.id() << std::endl;

    for (int j = 0; j < gd.n_polygons_with_holes(); j++) {
      std::cout << "Polygon " << j << " in GeoDiv" << std::endl;
      Polygon_with_holes pwh = gd.polygons_with_holes()[j];
      CGAL::Bbox_2 bb = pwh.bbox();
      for (double k = (unsigned int) bb.xmin() + 0.5; k < bb.xmax(); k++) {
        std::cout << "k = " << k << "\n";
        for (double l = (unsigned int) bb.ymin() + 0.5; l < bb.ymax(); l++) {

          //std::cout << "k = " << k << ", l = " << l << std::endl;

          // Test if coordinates (k, l) is in j-th polygon (possibly with
          // holes) of geographic division gd. If yes, then set density.
          if (CGAL::oriented_side(Point(k, l), pwh) == 1) {
            rho_init((int) k,(int) l) = 1.0;//map_state->target_areas_at(gd.id());
          }
        }
      }
    }
  }

  std::cout << "Breakpoint - Finished Setting Density" << std::endl;

  map_state->execute_fwd_plan();
  return;

}
