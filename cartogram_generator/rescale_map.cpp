#include <iostream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <list>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Cartesian.h>
#include <vector>
#include "constants.h"
#include "map_state.h"

typedef CGAL::Cartesian<double> K;
typedef CGAL::Aff_transformation_2<K> Transformation;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Vector_2<K> Vector;

void rescale_map(int long_lattice_side_length, MapState *map_state)
{

  double padding = (map_state->is_world_map() ?  1.0 : padding_unless_world);

  // Initialize bounding box of map with bounding box of 0-th
  // Polygon_with_holes in 0-th GeoDiv
  GeoDiv gd0 = map_state->get_geo_divs()[0];
  std::vector<Polygon_with_holes> pwhs = gd0.get_polygons_with_holes();
  CGAL::Bbox_2 bb0 = pwhs[0].bbox();
  double map_xmin = bb0.xmin();
  double map_xmax = bb0.xmax();
  double map_ymin = bb0.ymin();
  double map_ymax = bb0.ymax();

  // Expand bounding box to enclose all GeoDivs
  for (auto gd : map_state->get_geo_divs()) {
    for (auto pwh : gd.get_polygons_with_holes()) {
      CGAL::Bbox_2 bb = pwh.bbox();
      map_xmin = (bb.xmin() < map_xmin ? bb.xmin() : map_xmin);
      map_ymin = (bb.ymin() < map_ymin ? bb.ymin() : map_ymin);
      map_xmax = (bb.xmax() > map_xmax ? bb.xmax() : map_xmax);
      map_ymax = (bb.ymax() > map_ymax ? bb.ymax() : map_ymax);
    }
  }


  // Expand bounding box to guarantee a minimum padding
  double new_xmin = 0.5 * ((1.0-padding)*map_xmax + (1.0+padding)*map_xmin);
  double new_xmax = 0.5 * ((1.0+padding)*map_xmax + (1.0-padding)*map_xmin);
  double new_ymin = 0.5 * ((1.0-padding)*map_ymax + (1.0+padding)*map_ymin);
  double new_ymax = 0.5 * ((1.0+padding)*map_ymax + (1.0-padding)*map_ymin);

  // Ensure that the lattice dimensions lx and ly are integer powers of 2
  int lx, ly;
  double latt_const;
  if (map_xmax-map_xmin > map_ymax-map_ymin) {
    lx = long_lattice_side_length;
    latt_const = (new_xmax-new_xmin) / lx;
    ly = 1 << ((int) ceil(log2((new_ymax-new_ymin) / latt_const)));
    new_ymax = 0.5*(map_ymax+map_ymin) + 0.5*ly*latt_const;
    new_ymin = 0.5*(map_ymax+map_ymin) - 0.5*ly*latt_const;
  } else {
    ly = long_lattice_side_length;
    latt_const = (new_ymax-new_ymin) / ly;
    lx = 1 << ((int) ceil(log2((new_xmax-new_xmin) / latt_const)));
    new_xmax = 0.5*(map_xmax+map_xmin) + 0.5*lx*latt_const;
    new_xmin = 0.5*(map_xmax+map_xmin) - 0.5*lx*latt_const;
  }
  std::cerr << "Using a " << lx << "-by-" << ly
            << " lattice with bounding box" << std::endl;
  std::cerr << "\t("
            << new_xmin << ", " << new_ymin << ", "
            << new_xmax << ", " << new_ymax << ")"
            << std::endl;

  // Set lattice dimensions in map_state
  map_state->set_lx(lx);
  map_state->set_ly(ly);

  // Rescale all GeoDiv coordinates
  Transformation translate(CGAL::TRANSLATION,
                           CGAL::Vector_2<Epick>(-new_xmin, -new_ymin));
  Transformation scale(CGAL::SCALING, (1.0/latt_const));
  for (auto &gd : *map_state->ref_to_geo_divs()) {
    for (auto &pwh : *gd.ref_to_polygons_with_holes()) {
      Polygon *ext_ring = &pwh.outer_boundary();
      *ext_ring = transform(translate, *ext_ring);
      *ext_ring = transform(scale, *ext_ring);
      for (auto hi = pwh.holes_begin(); hi != pwh.holes_end(); ++hi) {
        *hi = transform(translate, *hi);
        *hi = transform(scale, *hi);
      }
    }
  }
  return;
  
}
