#include "densify.h"
#include "helper_functions.h"

std::list<Polyline> densify(std::list<Polyline> polyline_list) {
  std::list<Polyline> polyline_list_dens; 
  for (Polyline polyl : polyline_list) {
    Polyline polyl_dens = polyl;

    /*
    polyl_dens.push_back(polyl[0]); 
    for (int i = 0; i < (int) polyl.size() - 1; i++) {
      Point a = polyl[i];
      Point b = polyl[i + 1];

      std::vector<Point> dens_pts = graticule_intersections(a, b);
      for (Point pt : dens_pts)
        polyl_dens.push_back(pt);
      polyl_dens.push_back(b);
    }*/

    std::cout << polyl.size() << " " << polyl_dens.size() << std::endl;
    polyline_list_dens.push_back(polyl_dens);
  }
  return polyline_list_dens;
}

