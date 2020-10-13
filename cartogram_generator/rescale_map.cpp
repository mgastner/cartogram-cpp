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
  double new_xmin = 0.5 * ((1.0-padding)*map_xmax + (1.0+padding)*map_xmin);
  double new_ymin = 0.5 * ((1.0-padding)*map_ymax + (1.0+padding)*map_ymin);
  double new_xmax = 0.5 * ((1.0+padding)*map_xmax + (1.0-padding)*map_xmin);
  double new_ymax = 0.5 * ((1.0+padding)*map_ymax + (1.0-padding)*map_ymin);
  int lx, ly;
  if (map_xmax-map_xmin > map_ymax-map_ymin) {
    lx = longer_lattice_length;
    double latt_const = (new_xmax-new_xmin) / longer_lattice_length;
    ly = 1 << ((int) ceil(log2((new_ymax-new_ymin) / latt_const)));
    new_ymax = 0.5*(map_ymax+map_ymin) + 0.5*ly*latt_const;
    new_ymin = 0.5*(map_ymax+map_ymin) - 0.5*ly*latt_const;
  }

  map_state->set_lx(lx);
  map_state->set_ly(ly);
  return;
}
