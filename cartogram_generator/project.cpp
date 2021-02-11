#include <boost/multi_array.hpp>
#include <iostream>
#include <vector>
#include "interpol.h"

#include "project.h"

void project(MapState *map_state){
  
  const unsigned int lx = map_state->lx();
  const unsigned int ly = map_state->ly();
  
  boost::multi_array<XYPoint, 2> &proj = *map_state->proj();
  
  // Calculate displacement from proj array
  boost::multi_array<double, 2> xdisp(boost::extents[lx][ly]);
  boost::multi_array<double, 2> ydisp(boost::extents[lx][ly]);

  for (unsigned int i = 0; i < lx; i++){
    for (unsigned int j=0; j<ly; j++) {
      xdisp[i][j] = proj[i][j].x - i - 0.5;
      ydisp[i][j] = proj[i][j].y - j - 0.5;
    }
  }

  std::vector<GeoDiv> new_geo_divs;

  for (auto gd : map_state->geo_divs()){
    // For each GeoDiv
    GeoDiv new_gd(gd.id());

    for (auto pwh : gd.polygons_with_holes()){
      // For each polygon with holes

      Polygon old_ext_ring = pwh.outer_boundary();
      Polygon new_ext_ring;

      for (unsigned int i = 0; i < old_ext_ring.size(); i++){
        // Update exterior ring coordinates
        new_ext_ring.push_back(Point(
          interpol(old_ext_ring[i][0], old_ext_ring[i][1],
                   &xdisp, 'x', lx, ly) + old_ext_ring[i][0],
          interpol(old_ext_ring[i][0], old_ext_ring[i][1],
                   &ydisp, 'y', lx, ly) + old_ext_ring[i][1]
        ));
      }

      std::vector<Polygon> hole_v;

      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); hci++){
        Polygon old_hole = *hci;
        Polygon new_hole;

        for (unsigned int i = 0; i < old_hole.size(); i++){
          // Update hole coordinates
          new_hole.push_back(Point(
            interpol(old_hole[i][0], old_hole[i][1],
                     &xdisp, 'x', lx, ly) + old_hole[i][0],
            interpol(old_hole[i][0], old_hole[i][1],
                     &ydisp, 'y', lx, ly) + old_hole[i][1]
          ));
        }

        hole_v.push_back(new_hole);
      }

      const Polygon_with_holes new_pwh(new_ext_ring,
                                       hole_v.begin(),
                                       hole_v.end());
      new_gd.push_back(new_pwh);
    }
    new_geo_divs.push_back(new_gd);
  }
  map_state->set_geo_divs(new_geo_divs);
}