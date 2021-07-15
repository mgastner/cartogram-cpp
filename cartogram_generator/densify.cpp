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

    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      Polygon outer = pgnwh.outer_boundary();
      Polygon outer_dens;
      size_t outer_size = outer.size();

      // Iterate over each point in the outer boundary of the GeoDiv
      for (size_t i = 0; i < outer_size; i++) {
        
        // The segment defined by points a and b is to be densified.
        // b should be the point immediately after a, unless a is the final
        // point of the boundary, in which case b should be the first point.
        Point a = outer[i];
        Point b = (i == outer_size - 1) ? outer[0] : outer[i + 1];

        // Densify the segment.
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

        // Pushing all points except the last.
        for (size_t i = 0; i < (outer_pts_dens.size() - 1); i++) {
          outer_dens.push_back(outer_pts_dens[i]);
        }
      }

      // Check for duplicate points in the densified outer boundary
      std::vector<Point> temp_out;
      for (size_t i = 0; i < outer_dens.size(); i++) {
        temp_out.push_back(outer_dens[i]);
      }
      if (duplicates(temp_out)) {
        std::cout << "Duplicates found in outer boundary!" << std::endl;
      }

      std::vector<Polygon> holes_v_dens;
      std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());

      // Iterate over each hole
      for (Polygon hole : holes_v) {
        Polygon hole_dens;
        for (size_t j = 0; j < hole.size(); j++) {
          
          // c and d are determined the same way as a and b above
          Point c = hole[j];
          Point d = (j == hole.size() - 1) ? hole[0] : hole[j + 1];

          std::vector<Point> hole_pts_dens = densification_points(c, d);

          if (abs(c[0] - d[0]) > 10 || abs(c[1] - d[1]) > 10) {
            std::cout << "Points VERY different!" << std::endl;
            std::cout << "c: " << c << std::endl;
            std::cout << "d: " << d << std::endl;
          }

          // if (duplicates(hole_pts_dens)) {
          //   std::cout << "Original: c " << c << " d " << d << std::endl;
          //   std::cout << "Densified: " << std::endl;
          //   for (size_t i = 0; i < hole_pts_dens.size(); i++) {
          //     std::cout << i + 1 << ": ";
          //     std::cout << hole_pts_dens[i] << std::endl;
          //   }
          // }

          for (size_t i = 0; i < (hole_pts_dens.size() - 1); i++) {
            hole_dens.push_back(hole_pts_dens[i]);
          }
        }

        // Check for duplicate points in the densified hole boundary
        std::vector<Point> temp_holes;
        for (size_t i = 0; i < hole_dens.size(); i++) {
          temp_holes.push_back(hole_dens[i]);
        }
        if (duplicates(temp_holes)) {
          std::cout << "Duplicates found in hole!" << std::endl;
        }

        holes_v_dens.push_back(hole_dens);
      }
      Polygon_with_holes pgnwh_dens(outer_dens, holes_v_dens.begin(), holes_v_dens.end());
      gd_dens.push_back(pgnwh_dens);
    }
    geodivs_dens.push_back(gd_dens);
  }

  std::cout << "Finished densifying\n";
  return geodivs_dens;
}
