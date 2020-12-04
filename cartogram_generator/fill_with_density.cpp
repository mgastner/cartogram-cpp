#include "map_state.h"
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Bbox_2.h>

void fill_with_density(MapState * map_state) {
  double total_current_area = 0.0;
  for (auto gd: map_state -> geo_divs()) {
    total_current_area += gd.area();
  }
  double total_target_area = 0.0;
  for (auto gd: map_state -> geo_divs()) {
    total_target_area += map_state -> target_areas_at(gd.id());
  }

  double mean_density = total_target_area / total_current_area;
  FTReal2d & rho_init = * map_state -> ref_to_rho_init();

  for (unsigned int i = 0; i < map_state -> lx(); i++) {
    for (unsigned int j = 0; j < map_state -> ly(); j++) {
      rho_init(i, j) = 0; // set those that are 0 to mean density later
    }
  }

  // Resolution of density
  double res = 16;

  // Iterating through geodivs in map_state
  for (auto gd: map_state -> geo_divs()) {
    std::cout << "Working on gd with ID: " << gd.id() << std::endl;

    // Associative area only called once
    double target_density = map_state -> target_areas_at(gd.id()) / gd.area();
    target_density /= res;

    // Iterating through polygons with holes in map_state
    for (int j = 0; j < gd.n_polygons_with_holes(); j++) {
      std::cout << "Polygon " << j << " in GeoDiv" << std::endl;
      Polygon_with_holes pwh = gd.polygons_with_holes()[j];
      CGAL::Bbox_2 bb = pwh.bbox();

      // cylce through y co-ordinates
      for (double k = (unsigned int) floor(bb.ymin()) - 0.5; k < ceil(bb.ymax()) + 0.5; k++) {

        for (double line_y = k; line_y < k + 1; line_y += (1.0 / res)) {
          Polygon ext_ring = pwh.outer_boundary();
          double prev_point[2];
          prev_point[0] = ext_ring[0][0];
          prev_point[1] = ext_ring[0][1];
          std::vector < double > intersections;

          // Run algorithm on exterior ring
          for (size_t l = 1; l < ext_ring.size(); l++) {
            double curr_point[2];
            curr_point[0] = ext_ring[l][0];
            curr_point[1] = ext_ring[l][1];
            // Intersection present
            if ((curr_point[1] <= line_y && prev_point[1] >= line_y) ||
              (curr_point[1] >= line_y && prev_point[1] <= line_y)) {

              if (curr_point[1] == prev_point[1]) {
                std::cout << "Grazing Incident" << std::endl;
                curr_point[1] += (curr_point[1] * 1/res);
                prev_point[1] -= (prev_point[1] * 1/res);
              }

              if (curr_point[1] == line_y || prev_point[1] == line_y) {
                std::cout << '\n';
                std::cout << "THIS HAS RUN" << '\n';
                std::cout << '\n';
                l++; // otherwise, the same point will be counted twice
              }

              double intersection = (curr_point[0] * (prev_point[1] - line_y) +
                  prev_point[0] * (line_y - curr_point[1])) /
                (prev_point[1] - curr_point[1]);
              intersections.push_back(intersection);
            }
            prev_point[0] = curr_point[0];
            prev_point[1] = curr_point[1];
          }

          // Run algorithm on each hole
          for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
            Polygon hole = * hci;
            prev_point[0] = hole[0][0];
            prev_point[1] = hole[0][1];

            for (size_t l = 1; l < hole.size(); l++) {

              double curr_point[2];
              curr_point[0] = hole[l][0];
              curr_point[1] = hole[l][1];
              // intersection present
              if ((curr_point[1] <= line_y && prev_point[1] >= line_y) ||
                (curr_point[1] >= line_y && prev_point[1] <= line_y)) {

                if (curr_point[1] == prev_point[1]) {
                  std::cout << "Grazing Incident" << std::endl;
                  curr_point[1] += (curr_point[1] * 1/res);
                  prev_point[1] -= (prev_point[1] * 1/res);
                }

                if (curr_point[1] == line_y || prev_point[1] == line_y) {
                  std::cout << '\n';
                  std::cout << "THIS HAS RUN" << '\n';
                  std::cout << '\n';
                  l++; // otherwise, the same point will be counted twice
                }

                double intersection = (curr_point[0] * (prev_point[1] - line_y) +
                    prev_point[0] * (line_y - curr_point[1])) /
                  (prev_point[1] - curr_point[1]);
                intersections.push_back(intersection);
              }
              prev_point[0] = curr_point[0];
              prev_point[1] = curr_point[1];
            }
          }

          if (intersections.size() % 2 != 0) {
            std::cout << "Incorrect Topology" << std::endl;
            std::cout << "Number of intersections: " << intersections.size();
            std::cout << std::endl;
            continue; // Ignoring for now
          }

          // sorting vector in ascending order of intersection
          sort(intersections.begin(), intersections.end());

          for (size_t l = 0; l < intersections.size(); l += 2) {
            double left_x = intersections[l];
            double right_x = intersections[l + 1];
            for (size_t m = ceil(left_x); m < right_x; m++) {
              if (m == ceil(left_x)) {
                rho_init((int) m, (int) ceil(k)) += ( target_density *
                                                    (ceil(left_x) - left_x) );
              }
              else if (m == floor(right_x)) {
                rho_init((int) m, (int) ceil(k)) += ( target_density *
                                                    (right_x - floor(right_x)) );
              }
              else {
                  rho_init((int) m, (int) ceil(k)) += target_density;
              }
            }
          }

        }
      }
    }
  }

  for (unsigned int i = 0; i < map_state -> lx(); i++) {
    for (unsigned int j = 0; j < map_state -> ly(); j++) {
      if (rho_init(i, j) == 0) {
        rho_init(i, j) = mean_density;
      }
    }
  }

  map_state -> execute_fwd_plan();
  return;

}
