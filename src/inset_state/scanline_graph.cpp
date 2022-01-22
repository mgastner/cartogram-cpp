#include "../inset_state.h"
#include "../constants.h"

std::vector<std::vector<intersection> >
  InsetState::scanlines_parallel_to_axis(char grid_side, unsigned int res) const
{
  unsigned int grid_length = (grid_side == 'x') ? lx() : ly();
  int n_rays = grid_length * res;
  std::vector<std::vector<intersection> > scanlines(n_rays);

  // Iterate through GeoDivs in inset_state
  for (auto gd : geo_divs()) {

    // Find target density
    double target_density;
    target_density = target_areas_.at(gd.id()) / gd.area();

    // Iterate through "polygons with holes" in inset_state
    for (const auto &pwh : gd.polygons_with_holes()) {
      Bbox bb = pwh.bbox();

      double min_lim, max_lim;
      if (grid_side == 'x') {
        min_lim = bb.ymin(); max_lim = bb.ymax();
      } else {
        min_lim = bb.xmin(); max_lim = bb.xmax();
      }

      // Cycle through coordinates in bounding box of pwh
      for (unsigned int k = floor(min_lim) - 1;
           k <= ceil(max_lim) + 1;
           ++k) {

        // Cycle through each of the "test rays" between the graticule lines
        // y = k and y = k+1
        for (double ray = k + (1.0/res)/2;
             ray < k + 1;
             ray += (1.0/res)) {
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
            intersection temp = (grid_side == 'x') ? intersection(true) :
              intersection(false);
            if (temp.ray_intersects(curr_point,
                                    prev_point,
                                    ray,
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
              intersection temp = (grid_side == 'x') ? intersection(true) :
                intersection(false);
              if (temp.ray_intersects(curr_point,
                                      prev_point,
                                      ray,
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
            std::cerr << "Incorrect Topology." << std::endl;
            std::cerr << "Number of intersections: " << intersections.size();
            std::cerr << std::endl;
            std::cerr << grid_side << "-coordinate: " << ray << std::endl;
            std::cerr << "Intersection points: " << std::endl;
            for (unsigned int l = 0; l < intersections.size(); ++l) {
              std::cerr <<
                ((grid_side == 'x') ?
                 intersections[l].x() :
                 intersections[l].y())
                        << std::endl;
            }
            std::cerr << std::endl << std::endl;
            _Exit(932875);
          }
          std::sort(intersections.begin(), intersections.end());

          // Assign directions to intersections and add sorted vector of
          // intersections to vector scanlines
          int index = round((ray - 0.5/res) * res);
          for (unsigned int l = 0; l < intersections.size(); ++l) {
            intersections[l].direction = (l%2 == 0);
            scanlines[index].push_back(intersections[l]);
          }
        }
      }
    }
  }
  return scanlines;
}

// Creates continuity/adjacency graph using horizontal and vertical scans above
void InsetState::create_contiguity_graph(unsigned int res)
{

  // Getting the chosen graph.
  for (char graph : {'h', 'v'}) {
    std::vector<std::vector<intersection> > scan_graph;
    unsigned int max_k;
    if (graph == 'h') {
      scan_graph = scanlines_parallel_to_axis('x', res);
      max_k= ly();
    } else {
      scan_graph = scanlines_parallel_to_axis('y', res);
      max_k= lx();
    }

    // Iterating through scanline graph
    for (unsigned int k = 0; k < max_k; ++k) {

      // Cycle through each of the "res" number of rays in one cell
      for (double ray = k + (1.0/res)/2;
           ray < k + 1;
           ray += (1.0/res)) {

        // Intersections for one ray
        std::vector<intersection> intersections =
          scan_graph[static_cast<int>(round(((ray - (1.0/res)/2.0) * res)))];

        // Sort vector in ascending order of intersection
        sort(intersections.begin(), intersections.end());

        int size = intersections.size() - 1;

        // Fill GeoDivs by iterating through intersections
        for (int l = 1; l < size; l += 2) {
          double coord_1 = intersections[l].x();
          double coord_2 = intersections[l + 1].x();
          std::string gd_1 = intersections[l].geo_div_id;
          std::string gd_2 = intersections[l + 1].geo_div_id;

          if (gd_1 != gd_2 && coord_1 == coord_2) {
            for (auto &gd : geo_divs_) {
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

// Returns line segments highlighting intersection points using scan graphs
const std::vector<Segment> InsetState::intersecting_segments(unsigned int res)
  const
{
  std::vector<Segment> int_segments;
  for (char graph : {'h', 'v'}) {

    // Getting appropriate scan graph
    std::vector<std::vector<intersection> > scans;
    unsigned int max_k;
    if (graph == 'h') {
      scans = scanlines_parallel_to_axis('x', res);
      max_k= ly();
    } else {
      scans = scanlines_parallel_to_axis('y', res);
      max_k= lx();
    }

    // Iterating through scanline graph
    for (unsigned int k = 0; k < max_k; ++k) {

      // Cycle through each of the "res" number of rays in one cell
      for (double ray = k + (1.0/res)/2;
          ray < k + 1;
          ray += (1.0/res)) {

        // Intersections for one ray
        std::vector<intersection> intersections =
          scans[static_cast<int>(round(((ray - (1.0/res)/2.0) * res)))];

        // Sort vector in ascending order of intersection
        std::sort(intersections.begin(), intersections.end());
        int size = intersections.size() - 1;

        // Check for intersections, if both are entering/exiting
        for (int l = 0; l < size; ++l) {
          if (intersections[l].direction == intersections[l + 1].direction &&
              intersections[l].direction &&
              l + 2 <= size) {
            Segment temp;
            if (graph == 'h' &&
                intersections[l + 1].x() != intersections[l + 2].x()) {
              temp = Segment(
                              Point(intersections[l + 1].x(), ray),
                              Point(intersections[l + 2].x(), ray)
                              );
              int_segments.push_back(temp);
            } else if (intersections[l + 1].y() != intersections[l + 2].y()) {
              temp = Segment(
                              Point(ray, intersections[l + 1].y()),
                              Point(ray, intersections[l + 2].y())
                              );
              int_segments.push_back(temp);
            }
          }
        }
      }
    }
  }
  return int_segments;
}
