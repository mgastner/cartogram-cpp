#include "inset_state.h"

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
  int n_rays = static_cast<int>(lx_ * res);
  std::vector<std::vector<intersection> > vertical_adj(n_rays);

  // Creating vertical adjacency graph
  for (auto gd : geo_divs_) {

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
            vertical_adj[index].push_back(intersections[l]);
          }
        }
      }
    }
  }

  // Setting Vertical Adjacency graph for automatic coloring
  return vertical_adj;
}
