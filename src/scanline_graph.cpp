#include "inset_state.h"

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

      // Pre-condition to ignore grazing incidence (i.e., a line segment along
      // the polygon is exactly on the test ray)
      (a.y != b.y)) {
    if (a.y == ray_y) {
      a.y += epsilon;
    } else if (b.y == ray_y) {
      b.y += epsilon;
    }

    // Edit intersection passed by reference
    // coord stores the x coordinate.
    temp->coord = (a.x * (b.y - ray_y) + b.x * (ray_y - a.y)) / (b.y - a.y);
    temp->target_density = target_density;
    temp->direction = false;  // Temporary value
    return true;
  }
  return false;
}

const std::vector<std::vector<intersection> >
  InsetState::horizontal_scans(unsigned int res) const
{
  int n_rays = static_cast<int>(this->ly_ * res);
  std::vector<std::vector<intersection> > horizontal_scans(n_rays);

  // Iterate through GeoDivs in inset_state
  for (auto gd : this->geo_divs_) {

    // Find target density
    double target_density;
    target_density = this->target_areas_.at(gd.id()) / gd.area();

    // Iterate through "polygons with holes" in inset_state
    for (unsigned int j = 0; j < gd.n_polygons_with_holes(); ++j) {
      Polygon_with_holes pwh = gd.polygons_with_holes()[j];
      Bbox bb = pwh.bbox();

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
          double epsilon = 1e-6 / res;

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

          // Check if the number of intersections is odd
          if (intersections.size() % 2 != 0) {
            std::cerr << "Incorrect Topology" << std::endl;
            std::cerr << "Number of intersections: " << intersections.size();
            std::cerr << std::endl;
            std::cerr << "Y-coordinate: " << ray_y << std::endl;
            std::cerr << "Intersection points: " << std::endl;
            for (unsigned int l = 0; l < intersections.size(); ++l) {
              std::cerr << intersections[l].coord << std::endl;
            }
            std::cerr << std::endl << std::endl;
            _Exit(932875);
          }
          std::sort(intersections.begin(), intersections.end());

          // Add sorted vector of intersections to vector horizontal_scans
          for (unsigned int l = 0; l < intersections.size(); ++l) {
            intersections[l].direction = (l%2 == 0);
            int index = round((ray_y - 0.5/res) * res);
            horizontal_scans[index].push_back(intersections[l]);
          }
        }
      }
    }
  }
  return horizontal_scans;
}

bool ray_x_intersects(XYPoint a,
                      XYPoint b,
                      double ray_x,
                      intersection *temp,
                      double target_density,
                      double epsilon)
{
  // Check if intersection is present
  if (((a.x <= ray_x && b.x >= ray_x) ||
       (a.x >= ray_x && b.x <= ray_x)) &&

      // Pre-condition to ignore grazing incidence (i.e. a line segment along
      // the polygon is exactly on the test ray)
      (a.x != b.x)) {
    if (a.x == ray_x) {
      a.x += epsilon;
    } else if (b.x == ray_x) {
      b.x += epsilon;
    }

    // Edit intersection passed by reference
    // coord stores the y coordinate.
    temp->coord = (a.y * (b.x - ray_x) + b.y * (ray_x - a.x)) / (b.x - a.x);
    temp->target_density = target_density;
    temp->direction = false;  // Temporary value
    return true;
  }
  return false;
}

// Creates a vector of intersections between each GeoDiv and each scanline.
const std::vector<std::vector<intersection> >
  InsetState::vertical_scans(unsigned int res) const
{

  // A vector to store the vertical adjacency graph.
  // Inspired by code for fill_with_density.cpp
  int n_rays = static_cast<int>(this->lx_ * res);
  std::vector<std::vector<intersection> > vertical_scans(n_rays);

  // Creating vertical adjacency graph
  for (auto gd : this->geo_divs_) {

    // Iterate through "polygons with holes" in inset_state
    for (unsigned int j = 0; j < gd.n_polygons_with_holes(); ++j) {
      Polygon_with_holes pwh = gd.polygons_with_holes()[j];
      Bbox bb = pwh.bbox();

      // Cycle through x-coordinates in bounding box of pwh
      for (unsigned int k = floor(bb.xmin()) - 1;
           k <= ceil(bb.xmax()) + 1;
           ++k) {

        // Cycle through each of the "test rays" between the graticule lines
        // x = k and x = k+1
        for (double ray_x = k + (1.0/res)/2;
             ray_x < k + 1;
             ray_x += (1.0/res)) {
          Polygon ext_ring = pwh.outer_boundary();
          XYPoint prev_point;
          prev_point.x = ext_ring[ext_ring.size()-1][0];
          prev_point.y = ext_ring[ext_ring.size()-1][1];


          // Temporary vector of intersections for this particular ray
          std::vector<intersection> intersections;

          double epsilon = 1e-6 * (1.0/res);

          // First we run the algorithm on the exterior ring
          for (unsigned int l = 0; l < ext_ring.size(); ++l) {
            XYPoint curr_point;
            curr_point.x = ext_ring[l][0];
            curr_point.y = ext_ring[l][1];
            intersection temp;
            if (ray_x_intersects(curr_point,
                                 prev_point,
                                 ray_x,
                                 &temp,
                                 0,
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
              if (ray_x_intersects(curr_point,
                                   prev_point,
                                   ray_x,
                                   &temp,
                                   0,
                                   epsilon)) {

                temp.geo_div_id = gd.id();
                intersections.push_back(temp);
              }
              prev_point.x = curr_point.x;
              prev_point.y = curr_point.y;
            }
          }

          // Check if odd number of intersections, indicating self-intersection
          if (intersections.size() % 2 != 0) {
            std::cerr << "Incorrect Topology" << std::endl;
            std::cerr << "Number of intersections: " << intersections.size();
            std::cerr << std::endl;
            std::cerr << "X-coordinate: " << ray_x << std::endl;
            std::cerr << "Intersection points: " << std::endl;
            for (unsigned int l = 0; l < intersections.size(); ++l) {
              std::cerr << intersections[l].coord << std::endl;
            }
            std::cerr << std::endl << std::endl;
            _Exit(932875);
          }
          std::sort(intersections.begin(), intersections.end());

          // Add sorted vector of intersections to adjacency graph
          for (unsigned int l = 0; l < intersections.size(); ++l) {
            intersections[l].direction = (l%2 == 0);
            int index = round(((ray_x - (1.0/res)/2.0) * res));
            vertical_scans[index].push_back(intersections[l]);
          }
        }
      }
    }
  }
  return vertical_scans;
}

// Creates an adjacency graph using horizontal and vertical scans above
void InsetState::create_adjacency_graph(unsigned int res)
{

  // Getting the chosen graph.
  for (char graph : {'h', 'v'}) {
    std::vector<std::vector<intersection> > scan_graph;
    unsigned int max_k = 0;
    if (graph == 'h') {
      scan_graph = this->horizontal_scans(res);
      max_k = this->ly();
    } else if (graph == 'v') {
      scan_graph = this->vertical_scans(res);
      max_k = this->lx();
    }

    // Iterating through horizontal adjacency graph
    for (unsigned int k = 0; k < max_k; ++k) {

      // Cycle through each of the "res" number of rays in one cell
      for (double ray = k + (1.0/res)/2;
           ray < k + 1;
           ray += (1.0/res)) {

        // The intersections for one ray
        std::vector<intersection> intersections =
          scan_graph[static_cast<int>(round(((ray - (1.0/res)/2.0) * res)))];

        // Sort vector in ascending order of intersection
        sort(intersections.begin(), intersections.end());

        int size = intersections.size() - 1;

        // Fill GeoDivs by iterating through intersections
        for (int l = 1; l < size; l += 2) {
          double x_1 = intersections[l].coord;
          double x_2 = intersections[l + 1].coord;
          std::string gd_1 = intersections[l].geo_div_id;
          std::string gd_2 = intersections[l + 1].geo_div_id;

          if (gd_1 != gd_2 && x_1 == x_2) {
            for (auto &gd : this->geo_divs_) {
              if (gd.id() == gd_1) {
                gd.adjacent_to(gd_2);
              } else if (gd.id() == gd_2) {
                gd.adjacent_to(gd_1);
              }
            }
          }

        }
      }
    }
  }

}

// Returns line segments highlighting intersection points using scan graphs.
const std::vector<Segment> InsetState::intersections() const
{
  std::vector<Segment> intersections;
  // std::vector<Polygon_with_holes> all_pgwhs_in_inset;
  return intersections;
}
