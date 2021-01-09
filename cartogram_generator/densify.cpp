#include "densify.h"
#include "find_graticule_intersections.h"

std::vector<GeoDiv> densify(std::vector<GeoDiv> container) {

  std::vector<GeoDiv> container_dens; 
  for (GeoDiv gd : container) {
    GeoDiv gd_dens(gd.id());

    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      Polygon outer = pgnwh.outer_boundary();
      Polygon outer_dens;
      outer_dens.push_back(outer[0]);
      for (int i = 0; i < (int) outer.size() - 1; i++) {
        Point a = outer[i];
        Point b = outer[i + 1];

        std::vector<Point> outer_pts_dens = graticule_intersections(a, b);
        for (Point outer_pt : outer_pts_dens)
          outer_dens.push_back(outer_pt);

        outer_dens.push_back(b);
      }

      std::vector<Polygon> holes_v_dens;
      std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
      for (Polygon hole : holes_v) {
        Polygon hole_dens;
        hole_dens.push_back(hole[0]);
        for (int j = 0; j < (int) hole.size() - 1; j++) {

          Point c = hole[j];
          Point d = hole[j + 1];

          std::vector<Point> holes_pts_dens = graticule_intersections(c, d);
          for (Point hole_pt : holes_pts_dens)
            hole_dens.push_back(hole_pt);

          hole_dens.push_back(d);
        }
        holes_v_dens.push_back(hole_dens);
      }

      Polygon_with_holes pgnwh_dens(outer_dens, holes_v_dens.begin(), holes_v_dens.end());  

      gd_dens.push_back(pgnwh_dens);
    }

    //std::cout << polyl.size() << " " << polyl_dens.size() << std::endl;
    container_dens.push_back(gd_dens);
  }
  return container_dens;
}

