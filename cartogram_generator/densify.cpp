#include "densify.h"
#include "densification_points.h"

bool duplicates(std::vector<Point> v) {
  CGAL::set_pretty_mode(std::cout);
  for (size_t i = 0; i < v.size() - 1; i++) {
    if (point_almost_equal(v[i], v[i + 1])) {
      std::cout << "Point: " << i << ", v[i]: " << v[i] << std::endl;
      std::cout << "Point: " << i + 1 << ", v[i + 1]: " << v[i + 1] << std::endl;
      return true;
    }
  }
  return false;
}

std::vector<GeoDiv> densify(std::vector<GeoDiv> geodivs) {

  std::cout << "Started densifying" << std::endl;
  std::vector<GeoDiv> geodivs_dens;
  for (GeoDiv gd : geodivs) {
    GeoDiv gd_dens(gd.id());
    std::cout << "Densifying: " << gd.id() << std::endl;

    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      Polygon outer = pgnwh.outer_boundary();
      Polygon outer_dens;
      size_t outer_size = outer.size();
      for (size_t i = 0; i < outer_size; i++) {
        Point a = outer[i];
        Point b;

        if (i == outer_size - 1){
          b = outer[0];
        } else {
          b = outer[i + 1];
        }

        // std::cout << "Densifying outer pts \na: " << a << "\nb: " << b << "\ni: " << i  << "\nouter_size: " << outer_size << "\n\n";
        std::vector<Point> outer_pts_dens = densification_points(a, b);

        if (abs(a[0] - b[0]) > 10 || abs(a[1] - b[1]) > 10) {
          std::cout << "Points VERY different!" << std::endl;
          std::cout << "a: " << a << std::endl;
          std::cout << "b: " << b << std::endl;
        }

        // For printing duplicates upto an epsilon
        // if (duplicates(outer_pts_dens)) {
        //   std::cout << "Original: a " << a << " b " << b << std::endl;
        //   std::cout << "Densified: " << std::endl;
        //   for (size_t i = 0; i < outer_pts_dens.size(); i++) {
        //     std::cout << i + 1 << ": ";
        //     std::cout << outer_pts_dens[i] << std::endl;
        //   }
        // }

        // Pushing all points but not the last one
        for (size_t i = 0; i < (outer_pts_dens.size() - 1); i++) {
          outer_dens.push_back(outer_pts_dens[i]);
        }

      }

      std::vector<Point> temp_out;
      for (size_t i = 0; i < outer_dens.size(); i++) {
        temp_out.push_back(outer_dens[i]);
      }

      if (duplicates(temp_out)) {
        // for (size_t i = 0; i < temp_out.size(); i++) {
          std::cout << "Duplicates found in outer boundary!" << std::endl;
        // }
      }

      std::vector<Polygon> holes_v_dens;
      std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
      for (Polygon hole : holes_v) {
        Polygon hole_dens;
        for (size_t j = 0; j < hole.size(); j++) {
          Point c = hole[j];
          Point d;

          if (j == hole.size() - 1){
            d = hole[0];
          } else {
            d = hole[j + 1];
          }

          std::vector<Point> hole_pts_dens = densification_points(c, d);

          if (duplicates(hole_pts_dens)) {
            std::cout << "Original: c " << c << " d " << d << std::endl;
            std::cout << "Densified: " << std::endl;
            for (size_t i = 0; i < hole_pts_dens.size(); i++) {
              std::cout << i + 1 << ": ";
              std::cout << hole_pts_dens[i] << std::endl;
            }
          }

          for (size_t i = 0; i < (hole_pts_dens.size() - 1); i++) {
            hole_dens.push_back(hole_pts_dens[i]);
          }

        }

        std::vector<Point> temp_holes;
        for (size_t i = 0; i < hole_dens.size(); i++) {
          temp_holes.push_back(hole_dens[i]);
        }

        if (duplicates(temp_holes)) {
          // for (size_t i = 0; i < temp_out.size(); i++) {
            std::cout << "Duplicates found in hole!" << std::endl;
          // }
        }

        holes_v_dens.push_back(hole_dens);
      }

      Polygon_with_holes pgnwh_dens(outer_dens, holes_v_dens.begin(), holes_v_dens.end());

      gd_dens.push_back(pgnwh_dens);
    }

    std::cout << "Densified : " << gd.id() << std::endl;
    // std::cout << polyl.size() << " --> " << polyl_dens.size() << std::endl;
    geodivs_dens.push_back(gd_dens);
  }
  return geodivs_dens;
}
