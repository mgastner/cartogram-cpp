#include "cartogram_info.h"
#include "inset_state.h"
#include "auto_color.h"
#include "colors.h"

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
    // Although it's name is x, this is the y-coordinate. This is becasue
    // the intersections struct was initally made for fill_with_density.cpp
    // which only required us to run through horizontal rays.
    temp->x = (a.y * (b.x - ray_x) + b.y * (ray_x - a.x)) / (b.x - a.x);
    temp->target_density = target_density;
    temp->direction = false;  // Temporary value
    return true;
  }
  return false;
}

void create_vertical_adjacency_graph(InsetState* inset_state, unsigned int res)
{

  // A vector to store the vertical adjacency graph.
  // Inspired by code for fill_with_density.cpp
  int n_rays = (int) (inset_state->lx() * res);
  std::vector<std::vector<intersection> > vertical_adj(n_rays);

  // Creating vertical adjacency graph
  for (auto gd : inset_state->geo_divs()) {

    // Iterate through "polygons with holes" in inset_state
    for (int j = 0; j < gd.n_polygons_with_holes(); ++j) {
      Polygon_with_holes pwh = gd.polygons_with_holes()[j];
      CGAL::Bbox_2 bb = pwh.bbox();

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
              std::cerr << intersections[l].x << std::endl;
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
  inset_state->set_vertical_adj(vertical_adj);
}

void create_adjacency_graph(InsetState* inset_state,
                            char graph,
                            unsigned int res)
{

  // Getting the chosen graph.
  std::vector<std::vector<intersection> > adj_graph;
  unsigned int max_k = 0;
  if (graph == 'h') {
    adj_graph = inset_state->horizontal_adj();
    max_k = inset_state->ly();
  } else if (graph == 'v') {
    adj_graph = inset_state->vertical_adj();
    max_k = inset_state->lx();
  }

  // Iterating through horizontal adjacency graph
  for (unsigned int k = 0; k < max_k; ++k) {

    // Cycle through each of the "res" number of rays in one cell
    for (double ray = k + (1.0/res)/2;
         ray < k + 1;
         ray += (1.0/res)) {

      // The intersections for one ray
      std::vector<intersection> intersections =
        adj_graph[(int) round(((ray - (1.0/res)/2.0) * res))];

      // Sort vector in ascending order of intersection
      sort(intersections.begin(), intersections.end());

      int size = intersections.size() - 1;

      // Fill GeoDivs by iterating through intersections
      for (int l = 1; l < size; l += 2) {
        double x_1 = intersections[l].x;
        double x_2 = intersections[l + 1].x;
        std::string gd_1 = intersections[l].geo_div_id;
        std::string gd_2 = intersections[l + 1].geo_div_id;

        if (gd_1 != gd_2 && x_1 == x_2) {
          for (auto &gd : *inset_state->ref_to_geo_divs()) {
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

// Function to automatically color topology based on adjacency graph
void auto_color(InsetState* inset_state)
{

  std::vector<Color> palette;

  // Using default Palette for now
  // From https://colorbrewer2.org/
  palette.push_back(Color(27, 158, 119)); // turquoise
  palette.push_back(Color(117, 112, 179)); // purple
  palette.push_back(Color(231, 41, 138)); // pinkish-magenta
  palette.push_back(Color(102, 166, 30)); // green
  palette.push_back(Color(230, 171, 2)); // yellow
  palette.push_back(Color(217, 95, 2)); // redish-orange
  palette.push_back(Color(102, 102, 102)); // dark grey
  palette.push_back(Color(166, 118, 29)); // mustard-brown


  // Getting the horizontal adjacency graph.
  std::vector<std::vector<intersection> > horizontal_adj =
    inset_state->horizontal_adj();

  // Find resolution
  unsigned int res = horizontal_adj.size() / inset_state->ly();

  // Creating vertical adjacency graph
  create_vertical_adjacency_graph(inset_state, res);

  // Creating full adjacency graph based on vertical and horizontal graphs
  create_adjacency_graph(inset_state, 'h', res);
  create_adjacency_graph(inset_state, 'v', res);

  // Count to maximize colors used. This changes the starting color that
  // the algorithm choses for each GeoDiv.
  unsigned int count = 0;

  // In case coloring was unsucesfull, we restrict to a smaller count.
  unsigned int max_i = palette.size();

  // Iterating until we are able to color the entire map
  while (inset_state->colors_size() < inset_state->n_geo_divs()) {

    // Iterating through GeoDivs
    for (auto &gd : *inset_state->ref_to_geo_divs()) {

      // Iterating through all possible colors
      for (size_t i = (count % max_i); i < palette.size(); i++) {
        Color c = palette[i];
        bool shared_color = false;

        // Iterating through adjacent GeoDivs to check whether shared color
        for (std::string gd_id : gd.adjacent_geodivs()) {

          // Checking to see whether color found
          if (inset_state->color_found(gd_id)) {
            if (inset_state->colors_at(gd_id) == c) {
              shared_color = true;
            }
          }
        }

        // If the color is not shared with any other adjacent GeoDiv,
        // we assign it the color.
        if (!shared_color) {
          inset_state->colors_insert(gd.id(), c);
          i = palette.size();
        }
      }
      count++;
    }
    max_i--;
  }
  return;
}
