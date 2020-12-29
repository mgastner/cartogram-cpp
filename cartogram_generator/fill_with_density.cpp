#include "map_state.h"
#include "write_eps.h"

// Struct to store intersection data
struct intersection {
  double x;  // x-coordinate of intersection
  double target_density;  // GeoDiv's target_density
  bool direction;  // Does intersection enter or exit?

  // Overload "<" operator for this data type. Idea from
  // https://stackoverflow.com/questions/4892680/sorting-a-vector-of-structs
  bool operator < (const intersection &rhs) const
  {
    return (x < rhs.x || (x == rhs.x && direction < rhs.direction));
  }
};

void fill_with_density(MapState* map_state)
{
  double total_current_area = 0.0;
  for (auto gd : map_state->geo_divs()) {
    total_current_area += gd.area();
  }
  double total_target_area = 0.0;
  for (auto gd : map_state->geo_divs()) {
    total_target_area += map_state->target_areas_at(gd.id());
  }
  double mean_density = total_target_area / total_current_area;
  FTReal2d &rho_init = *map_state->ref_to_rho_init();

  // Initially assign 0 to all densities
  for (unsigned int i = 0; i < map_state->lx(); i++) {
    for (unsigned int j = 0; j < map_state->ly(); j++) {
      rho_init(i, j) = 0;
    }
  }

  // Resolution with which we sample polygons. "res" is the number of
  // horizontal "test lines" between each of the ly consecutive horizontal
  // graticule lines.
  unsigned int res = 16;

  // A vector (map_intersections) to store vectors of intersections
  int n_lines = (int) (map_state->ly() * res);
  std::vector<std::vector<intersection>> map_intersections(n_lines);

  // Iterate through GeoDivs in map_state
  for (auto gd : map_state->geo_divs()) {

    // Associative area. It is only called once to find out the target
    // density.
    double target_density = map_state->target_areas_at(gd.id()) / gd.area();
    target_density /= res;

    // Iterate through "polygons with holes" in map_state
    for (int j = 0; j < gd.n_polygons_with_holes(); j++) {
      Polygon_with_holes pwh = gd.polygons_with_holes()[j];
      CGAL::Bbox_2 bb = pwh.bbox();

      // Cycle through y-coordinates in bounding box of pwh
      for (unsigned int k = floor(bb.ymin()) - 1; k <= ceil(bb.ymax()) + 1;
           k++) {

        // Cycle through each of the "test lines" between the graticule lines
        // y = k and y = k+1.
        for (double line_y = k + (1.0/res)/2; line_y < k + 1;
             line_y += (1.0/res)) {
          Polygon ext_ring = pwh.outer_boundary();
          double prev_point[2];
          prev_point[0] = ext_ring[ext_ring.size()-1][0];
          prev_point[1] = ext_ring[ext_ring.size()-1][1];

          // Temporary vector of intersections for this particular line
          std::vector<intersection> intersections;

          // The following algorithm works by iterating through "res" lines in
          // each cell. For each line, we iterate through every edge in a
          // polygon and store any intersections. Finally, once all
          // intersections have been stored, we iterate between intersections,
          // and add the appropriate densities.
          // We add a small value "epsilon" in case the line with equation
          // y = line_y goes exactly through curr_point. The addition ensures
          // that, if there is any intersection, it is only counted once. It
          // also correctly detects whether the line crosses through the point
          // without entering or exiting the polygon.
          double epsilon = 1e-6;

          // First we run the algorithm on the exterior ring
          for (unsigned int l = 0; l < ext_ring.size(); l++) {
            double curr_point[2];
            curr_point[0] = ext_ring[l][0];
            curr_point[1] = ext_ring[l][1];

            // Check if intersection is present
            if (((curr_point[1] <= line_y && prev_point[1] >= line_y) ||
                 (curr_point[1] >= line_y && prev_point[1] <= line_y)) &&

                // Pre-condition to ignore grazing incidence (i.e. a line
                // segment along the polygon is exactly on the test line)
                (curr_point[1] != prev_point[1])) {
              if (curr_point[1] == line_y) {
                curr_point[1] += epsilon * (1.0/res);
              } else if (prev_point[1] == line_y) {
                prev_point[1] += epsilon * (1.0/res);
              }

              // Create an intersection and store it in a vector
              intersection temp;
              temp.x = (curr_point[0] * (prev_point[1] - line_y) +
                        prev_point[0] * (line_y - curr_point[1])) /
                       (prev_point[1] - curr_point[1]);
              temp.target_density = target_density;
              temp.direction = false;   // Temporary value
              intersections.push_back(temp);
            }
            prev_point[0] = curr_point[0];
            prev_point[1] = curr_point[1];
          }

          // Run algorithm on each hole
          for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
            Polygon hole = *hci;
            prev_point[0] = hole[hole.size()-1][0];
            prev_point[1] = hole[hole.size()-1][1];
            for (unsigned int l = 0; l < hole.size(); l++) {
              double curr_point[2];
              curr_point[0] = hole[l][0];
              curr_point[1] = hole[l][1];
              if (((curr_point[1] <= line_y && prev_point[1] >= line_y) ||
                   (curr_point[1] >= line_y && prev_point[1] <= line_y)) &&
                  (curr_point[1] != prev_point[1])) {
                if (curr_point[1] == line_y) {
                  curr_point[1] += epsilon * (1.0/res);
                } else if (prev_point[1] == line_y) {
                  prev_point[1] += epsilon * (1.0/res);
                }
                intersection temp;
                temp.x = (curr_point[0] * (prev_point[1] - line_y) +
                          prev_point[0] * (line_y - curr_point[1])) /
                         (prev_point[1] - curr_point[1]);
                temp.target_density = target_density;
                temp.direction = false;  // Temporary value
                intersections.push_back(temp);
              }
              prev_point[0] = curr_point[0];
              prev_point[1] = curr_point[1];
            }
          }

          // Check if odd number of intersections
          if (intersections.size() % 2 != 0) {
            std::cerr << "Incorrect Topology" << std::endl;
            std::cerr << "Number of intersections: " << intersections.size();
            std::cerr << std::endl;
            std::cerr << "Y-coordinate: " << line_y << std::endl;
            std::cerr << "Intersection points: " << std::endl;
            for (unsigned int l = 0; l < intersections.size(); l++) {
              std::cerr << intersections[l].x << std::endl;
            }
            _Exit(932875);
          }
          std::sort(intersections.begin(), intersections.end());

          // Add sorted vector of intersections to vector map_intersections
          for (unsigned int l = 0; l < intersections.size(); l++) {
            intersections[l].direction = (l%2 == 0);
            int index = round(((line_y - (1.0/res)/2.0) * res));
            map_intersections[index].push_back(intersections[l]);
          }
        }
      }
    }
  }

  // Cycle through y-coordinates in map_state
  for (unsigned int k = 0; k < map_state->ly(); k++) {

    // Cycle through each of the "res" number of lines in one cell
    for (double line_y = k + (1.0/res)/2;
         line_y < k + 1;
         line_y += (1.0/res)) {
      std::vector<intersection> intersections =
        map_intersections[(int) round(((line_y - (1.0/res)/2.0) * res))];

      // Sort vector in ascending order of intersection
      sort(intersections.begin(), intersections.end());

      // Fill lines that have no intersections with mean_density
      if (intersections.size() == 0) {
        for (unsigned int l = 0; l < map_state->lx(); l++) {
          rho_init(l, k) += mean_density/res;
        }
      } else {

        // Fill from first coordinate up to first GeoDiv
        for (unsigned int l = 1; l <= ceil(intersections[0].x); l++) {
          if (l == ceil(intersections[0].x)) {
            rho_init(l - 1, k) +=
              (mean_density/res) *
              (intersections[0].x - floor(intersections[0].x));
          } else {
            rho_init(l - 1, k) += mean_density/res;
          }
        }

        // Fill any empty spaces between GeoDivs
        for (unsigned int l = 1; l < intersections.size() - 1; l += 2) {
          double left_x = intersections[l].x;
          double right_x = intersections[l + 1].x;

          // Pre-condition to ensure different intersecting points
          if (left_x != right_x) {
            for (unsigned int m = ceil(left_x); m <= ceil(right_x); m++) {
              if (ceil(left_x) == ceil(right_x)) {
                rho_init(m - 1, k) +=
                  intersections[l].target_density * (right_x - left_x);
              } else if (m == ceil(left_x)) {
                rho_init(m - 1, k) +=
                  ((mean_density/res) * (ceil(left_x) - left_x));
              } else if (m == ceil(right_x)) {
                rho_init(m - 1, k) +=
                  ((mean_density/res) * (right_x - floor(right_x)));
              } else {
                rho_init(m - 1, k) += (mean_density/res);
              }
            }
          }
        }

        // Fill from last GeoDiv up to last coordinate
        for (unsigned int l = ceil(intersections.back().x);
             l <= map_state->lx();
             l++) {
          if (l == ceil(intersections.back().x)) {
            rho_init(l - 1, k) +=
              (mean_density/res) *
              (ceil(intersections.back().x) - intersections.back().x);
          } else {
            rho_init(l - 1, k) += mean_density/res;
          }
        }
      }

      // Fill GeoDivs by iterating through intersections
      for (unsigned int l = 0; l < intersections.size(); l += 2) {
        double left_x = intersections[l].x;
        double right_x = intersections[l + 1].x;

        // Check for intersection of polygons, holes and GeoDivs
        if (intersections[l].direction == intersections[l + 1].direction) {

          // Highlight where intersection is present
          std::cerr << "Invalid Geometry!" << std::endl;
          std::cerr << "Intersection of Polygons/Holes/Geodivs" << std::endl;
          std::cerr << "Y-coordinate: " << k << std::endl;
          std::cerr << "Left X-coordinate: " << left_x << std::endl;
          std::cerr << "Right X-coordinate: " << right_x << std::endl;
          _Exit(8026519);
        }

        // Fill each cell between intersections
        for (unsigned int m = ceil(left_x); m <= ceil(right_x); m++) {
          if (ceil(left_x) == ceil(right_x)) {
            rho_init(m - 1, k) +=
              intersections[l].target_density * (right_x - left_x);
          } else if (m == ceil(left_x)) {
            rho_init(m - 1, k) +=
              intersections[l].target_density * (ceil(left_x) - left_x);
          } else if (m == ceil(right_x)) {
            rho_init(m - 1, k) +=
              intersections[l].target_density * (right_x - floor(right_x));
          } else {
            rho_init(m - 1, k) += intersections[l].target_density;
          }
        }
      }
    }
  }
  if (map_state->trigger_write_density_to_eps()) {
    std::string file_name =
      std::string("unblurred_density_") +
      std::to_string(map_state->n_finished_integrations()) +
      ".eps";
    std::cout << "Writing " << file_name << std::endl;
    write_density_to_eps(file_name, rho_init.array(), map_state);
  }
  map_state->execute_fwd_plan();
  return;
}
