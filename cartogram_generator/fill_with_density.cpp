#include "map_state.h"
#include "write_eps.h"
#include "fill_with_density.h"

bool line_y_intersects(XYPoint a,
                       XYPoint b,
                       double line_y,
                       intersection *temp,
                       double target_density,
                       double epsilon)
{
  // Check if intersection is present
  if (((a.y <= line_y && b.y >= line_y) ||
       (a.y >= line_y && b.y <= line_y)) &&

      // Pre-condition to ignore grazing incidence (i.e. a line segment along
      // the polygon is exactly on the test line)
      (a.y != b.y)) {
    if (a.y == line_y) {
      a.y += epsilon;
    } else if (b.y == line_y) {
      b.y += epsilon;
    }

    // Edit intersection passed by reference
    temp->x = (a.x * (b.y - line_y) + b.x * (line_y - a.y)) / (b.y - a.y);
    temp->target_density = target_density;
    temp->direction = false;  // Temporary value
    return true;
  }
  return false;
}


void fill_with_density(MapState* map_state)
{

  std::map<std::string, double> gd_to_number;

  for (GeoDiv gd : map_state->geo_divs()) {
    double temp = 0.0;
    gd_to_number.insert(std::pair<std::string, double>(gd.id(), temp));
  }

  // Calculate the total current area and total target area, excluding any
  // missing values
  double total_current_area = 0.0;
  for (auto gd : map_state->geo_divs()) {
    if (!map_state->target_area_is_missing(gd.id())) {
      total_current_area += gd.area();
    }
  }
  double total_target_area = 0.0;
  for (auto gd : map_state->geo_divs()) {
    if (!map_state->target_area_is_missing(gd.id())) {
      total_target_area += map_state->target_areas_at(gd.id());
    }
  }
  double mean_density = total_target_area / total_current_area;
  FTReal2d &rho_init = *map_state->ref_to_rho_init();

  // Initially assign 0 to all densities
  for (unsigned int i = 0; i < map_state->lx(); ++i) {
    for (unsigned int j = 0; j < map_state->ly(); ++j) {
      rho_init(i, j) = 0;
    }
  }

  // Resolution with which we sample polygons. "res" is the number of
  // horizontal "test lines" between each of the ly consecutive horizontal
  // graticule lines.
  unsigned int res = 16;

  // A vector (map_intersections) to store vectors of intersections
  int n_lines = (int) (map_state->ly() * res);
  std::vector<std::vector<intersection> > map_intersections(n_lines);

  // Iterate through GeoDivs in map_state
  for (auto gd : map_state->geo_divs()) {

    // Associative area. It is only called once to find out the target
    // density.
    double target_density;
    if (!map_state->target_area_is_missing(gd.id())) {
      target_density = map_state->target_areas_at(gd.id()) / gd.area();
      target_density /= res;
    } else {
      target_density = mean_density / res;
    }

    // Iterate through "polygons with holes" in map_state
    for (int j = 0; j < gd.n_polygons_with_holes(); ++j) {
      Polygon_with_holes pwh = gd.polygons_with_holes()[j];
      CGAL::Bbox_2 bb = pwh.bbox();

      // Cycle through y-coordinates in bounding box of pwh
      for (unsigned int k = floor(bb.ymin()) - 1;
           k <= ceil(bb.ymax()) + 1;
           ++k) {

        // Cycle through each of the "test lines" between the graticule lines
        // y = k and y = k+1
        for (double line_y = k + (1.0/res)/2;
             line_y < k + 1;
             line_y += (1.0/res)) {
          Polygon ext_ring = pwh.outer_boundary();
          XYPoint prev_point;
          prev_point.x = ext_ring[ext_ring.size()-1][0];
          prev_point.y = ext_ring[ext_ring.size()-1][1];


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
          double epsilon = 1e-6 * (1.0/res);

          // First we run the algorithm on the exterior ring
          for (unsigned int l = 0; l < ext_ring.size(); ++l) {
            XYPoint curr_point;
            curr_point.x = ext_ring[l][0];
            curr_point.y = ext_ring[l][1];
            intersection temp;
            if (line_y_intersects(curr_point,
                                  prev_point,
                                  line_y,
                                  &temp,
                                  target_density,
                                  epsilon)) {
              temp.geo_div_id = gd.id();
              intersections.push_back(temp);
            }
            prev_point.x = curr_point.x;
            prev_point.y = curr_point.y;
          }

          // Run algorithm on each hole
          for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
            Polygon hole = *hci;
            prev_point.x = hole[hole.size()-1][0];
            prev_point.y = hole[hole.size()-1][1];
            for (unsigned int l = 0; l < hole.size(); ++l) {
              XYPoint curr_point;
              curr_point.x = hole[l][0];
              curr_point.y = hole[l][1];
              intersection temp;
              if (line_y_intersects(curr_point,
                                    prev_point,
                                    line_y,
                                    &temp,
                                    target_density,
                                    epsilon)) {

                temp.geo_div_id = gd.id();
                intersections.push_back(temp);
              }
              prev_point.x = curr_point.x;
              prev_point.y = curr_point.y;
            }
          }

          // Check if odd number of intersections
          if (intersections.size() % 2 != 0) {
            std::cerr << "Incorrect Topology" << std::endl;
            std::cerr << "Number of intersections: " << intersections.size();
            std::cerr << std::endl;
            std::cerr << "Y-coordinate: " << line_y << std::endl;
            std::cerr << "Intersection points: " << std::endl;
            for (unsigned int l = 0; l < intersections.size(); ++l) {
              std::cerr << intersections[l].x << std::endl;
            }
            _Exit(932875);
          }
          std::sort(intersections.begin(), intersections.end());

          // Add sorted vector of intersections to vector map_intersections
          for (unsigned int l = 0; l < intersections.size(); ++l) {
            intersections[l].direction = (l%2 == 0);
            int index = round(((line_y - (1.0/res)/2.0) * res));
            map_intersections[index].push_back(intersections[l]);
          }
        }
      }
    }
  }

  // Cycle through y-coordinates in map_state
  for (unsigned int k = 0; k < map_state->ly(); ++k) {

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
        for (unsigned int l = 0; l < map_state->lx(); ++l) {
          rho_init(l, k) += mean_density/res;
        }
      } else {

        // Fill from first coordinate up to first GeoDiv
        for (unsigned int l = 1; l <= ceil(intersections[0].x); ++l) {
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
            for (unsigned int m = ceil(left_x); m <= ceil(right_x); ++m) {
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
             ++l) {
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
          // _Exit(8026519);
        }

        // Fill each cell between intersections
        for (unsigned int m = ceil(left_x); m <= ceil(right_x); ++m) {

          double td;

          std::string gd_id = intersections[l].geo_div_id;

          if (ceil(left_x) == ceil(right_x)) {
            rho_init(m - 1, k) +=
              intersections[l].target_density * (right_x - left_x);
            td = intersections[l].target_density * (right_x - left_x);
          } else if (m == ceil(left_x)) {
            rho_init(m - 1, k) +=
              intersections[l].target_density * (ceil(left_x) - left_x);
            td = intersections[l].target_density * (ceil(left_x) - left_x);
          } else if (m == ceil(right_x)) {
            rho_init(m - 1, k) +=
              intersections[l].target_density * (right_x - floor(right_x));
            td = intersections[l].target_density * (right_x - floor(right_x));
          } else {
            rho_init(m - 1, k) += intersections[l].target_density;
            td = intersections[l].target_density;
          }

          gd_to_number.at(gd_id) = gd_to_number.at(gd_id) + td;

        }
      }
    }
  }

  for (GeoDiv gd : map_state->geo_divs()) {
    std::cout << "ID: " << gd.id() << ", ";
    std::cout << "effective target area: "
              << gd_to_number.at(gd.id()) << '\n';
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
