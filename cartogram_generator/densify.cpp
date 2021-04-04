#include "densify.h"
#include "densification_points.h"

bool duplicates(std::vector<Point> v) {
  CGAL::set_pretty_mode(std::cout);
  for (size_t i = 0; i < v.size() - 1; i++) {
    if (point_almost_equal(v[i], v[i + 1])) {
      std::cout << "Point: " << i << ", v[i]: " << v[i] << '\n';
      std::cout << "Point: " << i + 1 << ", v[i + 1]: " << v[i + 1] << '\n';
      return true;
    }
  }
  return false;
}

std::vector<GeoDiv> densify(std::vector<GeoDiv> geodivs) {

  std::vector<GeoDiv> geodivs_dens;
  for (GeoDiv gd : geodivs) {
    GeoDiv gd_dens(gd.id());

    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      Polygon outer = pgnwh.outer_boundary();
      Polygon outer_dens;
      for (int i = 0; i < (int) outer.size(); i++) {
        // Point a = outer[i];
        // Point b = outer[i + 1];
        Point a;
        Point b;
        
        if (i == (int) outer.size() - 1){
          a = outer[i];
          b = outer[0];
        } else {
          a = outer[i];
          b = outer[i + 1];
        }

        std::vector<Point> outer_pts_dens = densification_points(a, b);

        if (duplicates(outer_pts_dens)) {
          std::cout << "Original: a " << a << " b " << b << '\n';
          std::cout << "Densified: " << '\n';
          for (size_t i = 0; i < outer_pts_dens.size(); i++) {
            std::cout << i + 1 << ": ";
            std::cout << outer_pts_dens[i] << '\n';
          }
        }

        for (size_t i = 0; i < (outer_pts_dens.size() - 1); i++) {
          outer_dens.push_back(outer_pts_dens[i]);
        }

        // if (i == (int) outer.size() - 2){
        //   for (size_t i = 0; i < outer_pts_dens.size(); i++) {
        //     outer_dens.push_back(outer_pts_dens[i]);
        //   }
        // } else {
        //   for (size_t i = 0; i < (outer_pts_dens.size() - 1); i++) {
        //     outer_dens.push_back(outer_pts_dens[i]);
        //   }
        // }

      }

      std::vector<Point> temp_out;
      for (size_t i = 0; i < outer_dens.size(); i++) {
        temp_out.push_back(outer_dens[i]);
      }

      if (duplicates(temp_out)) {
        // for (size_t i = 0; i < temp_out.size(); i++) {
          std::cout << "Duplicates found in outer boundary!" << '\n';
        // }
      }

      std::vector<Polygon> holes_v_dens;
      std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
      for (Polygon hole : holes_v) {
        Polygon hole_dens;
        for (int j = 0; j < (int) hole.size(); j++) {

          // Point c = hole[j];
          // Point d = hole[j + 1];

          Point c;
          Point d;

          if (j == (int) hole.size() - 1){
            c = hole[j];
            d = hole[0];
          } else {
            c = hole[j];
            d = hole[j + 1];
          }

          std::vector<Point> hole_pts_dens = densification_points(c, d);

          if (duplicates(hole_pts_dens)) {
            std::cout << "Original: c " << c << " d " << d << '\n';
            std::cout << "Densified: " << '\n';
            for (size_t i = 0; i < hole_pts_dens.size(); i++) {
              std::cout << i + 1 << ": ";
              std::cout << hole_pts_dens[i] << '\n';
            }
          }

          for (size_t i = 0; i < (hole_pts_dens.size() - 1); i++) {
            hole_dens.push_back(hole_pts_dens[i]);
          }

          // if (j == (int) outer.size() - 2){
          //   for (size_t i = 0; i < hole_pts_dens.size(); i++) {
          //     hole_dens.push_back(hole_pts_dens[i]);
          //   }
          // } else {
          //   for (size_t i = 0; i < (hole_pts_dens.size() - 1); i++) {
          //     hole_dens.push_back(hole_pts_dens[i]);
          //   }
          // }

        }

        std::vector<Point> temp_holes;
        for (size_t i = 0; i < hole_dens.size(); i++) {
          temp_holes.push_back(hole_dens[i]);
        }

        if (duplicates(temp_holes)) {
          // for (size_t i = 0; i < temp_out.size(); i++) {
            std::cout << "Duplicates found in hole!" << '\n';
          // }
        }

        holes_v_dens.push_back(hole_dens);
      }

      Polygon_with_holes pgnwh_dens(outer_dens, holes_v_dens.begin(), holes_v_dens.end());

      gd_dens.push_back(pgnwh_dens);
    }

    // std::cout << polyl.size() << " " << polyl_dens.size() << std::endl;
    geodivs_dens.push_back(gd_dens);
  }
  return geodivs_dens;
}
