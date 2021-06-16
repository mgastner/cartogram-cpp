#include "cartogram_info.h"
#include "inset_state.h"
#include "write_eps.h"
#include "fill_with_density.h"

bool ray_y_intersects(XYPoint a,
                      XYPoint b,
                      double ray_y,
                      intersection *temp,
                      double target_density,
                      double epsilon)
{
  // Check if intersection is present
  if (((a.y <= ray_y && b.y >= ray_y) ||
       (a.y >= ray_y && b.y <= ray_y)) &&

      // Pre-condition to ignore grazing incidence (i.e. a line segment along
      // the polygon is exactly on the test ray)
      (a.y != b.y)) {
    if (a.y == ray_y) {
      a.y += epsilon;
    } else if (b.y == ray_y) {
      b.y += epsilon;
    }

    // Edit intersection passed by reference
    temp->x = (a.x * (b.y - ray_y) + b.x * (ray_y - a.y)) / (b.y - a.y);
    temp->target_density = target_density;
    temp->direction = false;  // Temporary value
    return true;
  }
  return false;
}

void fill_with_density(InsetState* inset_state,
                       bool trigger_write_density_to_eps)
{

  // Calculate the total current area and total target area, excluding any
  // missing values
  double total_current_area = 0.0;
  for (auto gd : inset_state->geo_divs()) {
    if (!inset_state->target_area_is_missing(gd.id())) {
      total_current_area += gd.area();
    }
  }
  double total_target_area = 0.0;
  for (auto gd : inset_state->geo_divs()) {
    if (!inset_state->target_area_is_missing(gd.id())) {
      total_target_area += inset_state->target_areas_at(gd.id());
    }
  }
  double mean_density = total_target_area / total_current_area;
  FTReal2d &rho_init = *inset_state->ref_to_rho_init();

  // Initially assign 0 to all densities
  for (unsigned int i = 0; i < inset_state->lx(); ++i) {
    for (unsigned int j = 0; j < inset_state->ly(); ++j) {
      rho_init(i, j) = 0;
    }
  }

  // Resolution with which we sample polygons. "res" is the number of
  // horizontal "test rays" between each of the ly consecutive horizontal
  // graticule lines.
  unsigned int res = 16;

  // A vector (map_intersections) to store vectors of intersections
  int n_rays = (int) (inset_state->ly() * res);
  std::vector<std::vector<intersection> > map_intersections(n_rays);

  // Density numerator and denominator for each graticule cell
  // A density of a graticule cell can be calculated with (rho_num / rho_den)
  // Initialising all elements to 0 because we initially assume all graticule
  // cells to not be inside any GeoDiv. Any graticule cell where rho_den is 0
  // will get the mean_density
  std::vector<std::vector<double> >
    rho_num(inset_state->lx(), std::vector<double> (inset_state->ly(), 0));
  std::vector<std::vector<double> >
    rho_den(inset_state->lx(), std::vector<double> (inset_state->ly(), 0));

  // Iterate through GeoDivs in inset_state
  for (auto gd : inset_state->geo_divs()) {

    // Associative area. It is only called once to find out the target
    // density.
    double target_density;
    if (!inset_state->target_area_is_missing(gd.id())) {
      target_density = inset_state->target_areas_at(gd.id()) / gd.area();
    } else {
      target_density = mean_density;
    }

    // Iterate through "polygons with holes" in inset_state
    for (int j = 0; j < gd.n_polygons_with_holes(); ++j) {
      Polygon_with_holes pwh = gd.polygons_with_holes()[j];
      CGAL::Bbox_2 bb = pwh.bbox();

      // Cycle through y-coordinates in bounding box of pwh
      for (unsigned int k = floor(bb.ymin()) - 1;
           k <= ceil(bb.ymax()) + 1;
           ++k) {

        // Cycle through each of the "test rays" between the graticule lines
        // y = k and y = k+1
        for (double ray_y = k + (1.0/res)/2;
             ray_y < k + 1;
             ray_y += (1.0/res)) {
          Polygon ext_ring = pwh.outer_boundary();
          XYPoint prev_point;
          prev_point.x = ext_ring[ext_ring.size()-1][0];
          prev_point.y = ext_ring[ext_ring.size()-1][1];


          // Temporary vector of intersections for this particular ray
          std::vector<intersection> intersections;

          // The following algorithm works by iterating through "res" rays in
          // each cell. For each ray, we iterate through every edge in a
          // polygon and store any intersections. Finally, once all
          // intersections have been stored, we iterate between intersections,
          // and add the appropriate densities.
          // We add a small value "epsilon" in case the ray with equation
          // y = ray_y goes exactly through curr_point. The addition ensures
          // that, if there is any intersection, it is only counted once. It
          // also correctly detects whether the ray crosses through the point
          // without entering or exiting the polygon.
          double epsilon = 1e-6 * (1.0/res);

          // First we run the algorithm on the exterior ring
          for (unsigned int l = 0; l < ext_ring.size(); ++l) {
            XYPoint curr_point;
            curr_point.x = ext_ring[l][0];
            curr_point.y = ext_ring[l][1];
            intersection temp;
            if (ray_y_intersects(curr_point,
                                 prev_point,
                                 ray_y,
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
              if (ray_y_intersects(curr_point,
                                   prev_point,
                                   ray_y,
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
            std::cerr << "Y-coordinate: " << ray_y << std::endl;
            std::cerr << "Intersection points: " << std::endl;
            for (unsigned int l = 0; l < intersections.size(); ++l) {
              std::cerr << intersections[l].x << std::endl;
            }
            std::cerr << std::endl << std::endl;
            _Exit(932875);
          }
          std::sort(intersections.begin(), intersections.end());

          // Add sorted vector of intersections to vector map_intersections
          for (unsigned int l = 0; l < intersections.size(); ++l) {
            intersections[l].direction = (l%2 == 0);
            int index = round(((ray_y - (1.0/res)/2.0) * res));
            map_intersections[index].push_back(intersections[l]);
          }
        }
      }
    }
  }

  // Setting Horizontal Adjacency graph for automatic coloring
  inset_state->set_horizontal_adj(map_intersections);

  // Filling rho_num and rho_den (rho numerator and denominator)
  // rho_num is the sum of the weight * target_density for each segment of a
  // ray that is inside a GeoDiv
  // rho_den is the sum of the weights of a ray that is inside a GeoDiv
  // The weight of a segment of a ray that is inside a GeoDiv is calculated by
  // (length of the segment inside the geo_div) *
  // (area error of the geodiv)
  // Cycle through y-coordinates in inset_state
  for (unsigned int k = 0; k < inset_state->ly(); ++k) {

    // Cycle through each of the "res" number of rays in one cell
    for (double ray_y = k + (1.0/res)/2;
         ray_y < k + 1;
         ray_y += (1.0/res)) {

      // The intersections for one ray
      std::vector<intersection> intersections =
        map_intersections[(int) round(((ray_y - (1.0/res)/2.0) * res))];

      // Sort vector in ascending order of intersection
      sort(intersections.begin(), intersections.end());

      // Ensuring that the ray actually has intersections
      if (intersections.size() != 0) {

        // Fill any empty spaces between GeoDivs
        // where at least some part of a segment is inside a GeoDiv
        // (when a ray intersects a GeoDiv)
        for (unsigned int l = 1; l < intersections.size() - 1; l += 2) {
          double left_x = intersections[l].x;
          double right_x = intersections[l + 1].x;

          // Pre-condition to ensure different intersecting points
          if (left_x != right_x) {

            // We are intersecting with a GeoDiv. Hence, part of the ray
            // inside the graticule cell is inside the GeoDiv
            if (ceil(left_x) == ceil(right_x)) {
              double weight =
                inset_state->area_errs_at(intersections[l].geo_div_id) *
                                                        (right_x - left_x);
              double target_dens = intersections[l].target_density;
              rho_num[ceil(left_x) - 1][k] += weight * target_dens;
              rho_den[ceil(left_x) - 1][k] += weight;
            }
          }
        }

        // Fill last exiting intersection with GeoDiv where part of ray inside
        // the graticule cell is inside the GeoDiv
        unsigned int last_x = intersections.back().x;
        double last_weight =
          inset_state->area_errs_at(intersections.back().geo_div_id) *
                    (ceil(last_x) - last_x);
        double last_target_density = intersections.back().target_density;
        rho_num[ceil(last_x) - 1][k] += last_weight * last_target_density;
        rho_den[ceil(last_x) - 1][k] += last_weight;
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
          std::cerr << "Y-coordinate: " << ray_y << std::endl;
          std::cerr << "Left X-coordinate: " << left_x << std::endl;
          std::cerr << "Right X-coordinate: " << right_x << std::endl;
          std::cerr << std::endl;
          // _Exit(8026519);
        }

        // Fill each cell between intersections
        for (unsigned int m = ceil(left_x); m <= ceil(right_x); ++m) {

          double weight =
            inset_state->area_errs_at(intersections.back().geo_div_id);
          double target_dens = intersections[l].target_density;
          if (ceil(left_x) == ceil(right_x)) {
            weight *= (right_x - left_x);
          } else if (m == ceil(left_x)) {
            weight *= (ceil(left_x) - left_x);
          } else if (m == ceil(right_x)) {
            weight *= (right_x - floor(right_x));
          }
          rho_num[m - 1][k] += weight * target_dens;
          rho_den[m - 1][k] += weight;

        }
      }
    }
  }

  // Filling rho_init by dividing rho_num with rho_den
  for (unsigned int i = 0; i < inset_state->lx(); ++i) {
    for (unsigned int j = 0; j < inset_state->ly(); ++j) {
      if (rho_den[i][j] == 0) {
        rho_init(i, j) = mean_density;
      } else {
        rho_init(i, j) = rho_num[i][j] / rho_den[i][j];
      }
    }
  }

  if (trigger_write_density_to_eps) {
    std::string file_name =
      inset_state->inset_name() +
      "_unblurred_density_" +
      std::to_string(inset_state->n_finished_integrations()) +
      ".eps";
    std::cout << "Writing " << file_name << std::endl;
    write_density_to_eps(file_name, rho_init.array(), inset_state);
  }
  inset_state->execute_fwd_plan();
  return;
}
