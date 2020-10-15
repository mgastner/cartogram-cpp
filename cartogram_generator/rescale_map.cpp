#include <iostream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <list>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Cartesian.h>
#include <vector>

typedef CGAL::Cartesian<double> K;
typedef CGAL::Aff_transformation_2<K> Transformation;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Vector_2<K> Vector;

using namespace std;

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
  // declaring max values and change to be done with padding

  double latt_const, new_xmax, new_ymax, new_xmin, new_ymin;
  double lx, ly;
  double l = 8;

  /* Minimum dimensions that leave enough space between map and rectangular  */
  /* boundaries.                                                             */

  // https://github.com/mgastner/cartogram/blob/master/cartogram_generator/fill_with_density.c

  new_xmax = 0.5 * ((1.0+padding)*map_xmax + (1.0-padding)*map_xmin);
  new_xmin = 0.5 * ((1.0-padding)*map_xmax + (1.0+padding)*map_xmin);
  new_ymax = 0.5 * ((1.0+padding)*map_ymax + (1.0-padding)*map_ymin);
  new_ymin = 0.5 * ((1.0-padding)*map_ymax + (1.0+padding)*map_ymin);

  if (map_xmax-map_xmin > map_ymax-map_ymin) {
    lx = l;
    latt_const = (new_xmax-new_xmin) / l;
    cout << "Ly (before adjusting to power of 2): " << ly << endl;
    ly = 1 << ((int)ceil(log2((new_ymax-new_ymin)/latt_const)));
    cout << "Ly (after adjusting to power of 2): " << ly << endl;
    new_ymax = 0.5*(map_ymax+map_ymin) + 0.5*ly*latt_const;
    new_ymin = 0.5*(map_ymax+map_ymin) - 0.5*ly*latt_const;
  }
  else {
    ly = l;
    latt_const = (new_ymax-new_ymin) / l;
    lx = 1 << ((int) ceil(log2((new_xmax-new_xmin) / latt_const)));
    new_xmax = 0.5*(map_xmax+map_xmin) + 0.5*lx*latt_const;
    new_xmin = 0.5*(map_xmax+map_xmin) - 0.5*lx*latt_const;
  }

  cout << "Using a " << lx << " x " << ly  << " lattice with bounding box\n\t";
  cout << "(" << new_xmin << " " <<  new_ymin << " " << new_xmax << " ";
  cout << new_ymax << ")\n\n";

  /********************* Rescale all polygon coordinates. ********************/

  // PREVIOUS CODE
  //
  // for (i=0; i<n_poly; i++)
  //   for (j=0; j<n_polycorn[i]; j++) {
  //     polycorn[i][j].x = (polycorn[i][j].x - new_xmin) / latt_const;
  //     polycorn[i][j].y = (polycorn[i][j].y - new_ymin) / latt_const;
  //   }

  // creating transformation based on new values
  Transformation translate(CGAL::TRANSLATION, Vector(-new_xmin, -new_ymin));
  Transformation scale(CGAL::SCALING, (1/latt_const));

  // translating all points in polygon / transform every polygon
  // for (int i = 0; i < pgn.size(); i++) {
  //   cout << "original point: " << pgn[i] << "\n";
  //   pgn[i] = translate(pgn[i]);
  //   pgn[i] = scale(pgn[i]);
  //   cout << "translated point: " << pgn[i] << "\n\n";
  // }

  // // printing transformed polygon
  // cout << "translated polygon, pgn: \n " << pgn << endl;
  // cout << endl;

}
