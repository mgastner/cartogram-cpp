#include "constants.h"
#include "map_state.h"
#include <iostream>

void rescale_map(int longer_lattice_length, MapState *map_state)
{
  double padding = (map_state->is_world_map() ?  1.0 : padding_unless_world);

  // Initialize bounding box of map with bounding box of 0-th PolygonWH in
  // 0-th GeoDiv.
  GeoDiv gd0 = map_state->get_geo_divs()[0];
  std::vector<PolygonWH> pwhs = gd0.get_polygons_with_holes();
  CGAL::Bbox_2 bb0 = pwhs[0].bbox();
  double map_xmin = bb0.xmin();
  double map_xmax = bb0.xmax();
  double map_ymin = bb0.ymin();
  double map_ymax = bb0.ymax();
  for (auto gd : map_state->get_geo_divs()) {
    for (auto pwh : gd.get_polygons_with_holes()) {
      CGAL::Bbox_2 bb = pwh.bbox();
      map_xmin = (bb.xmin() < map_xmin ? bb.xmin() : map_xmin);
      map_xmax = (bb.xmax() < map_xmax ? bb.xmax() : map_xmax);
      map_ymin = (bb.ymin() < map_ymin ? bb.ymin() : map_ymin);
      map_ymax = (bb.ymax() < map_ymax ? bb.ymax() : map_ymax);
    }
  }
  std::cout << "Bbox of map: " << map_xmin << " " << map_ymin << " "
            << map_xmax << " " << map_ymax << std::endl;

  return;
}
