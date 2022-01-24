#include "../inset_state.h"
#include "../constants.h"

// This function adds intersections between a ray and a polygon to
// "intersections"
void add_intersections(std::vector<intersection> &intersections,
                       Polygon pgn,
                       double ray,
                       double target_density,
                       double epsilon,
                       std::string gd_id,
                       char grid_side)
{

          XYPoint prev_point;
          prev_point.x = pgn[pgn.size()-1][0];
          prev_point.y = pgn[pgn.size()-1][1];
          for (unsigned int l = 0; l < pgn.size(); ++l) {
            XYPoint curr_point;
            curr_point.x = pgn[l][0];
            curr_point.y = pgn[l][1];
            intersection temp = (grid_side == 'x') ? intersection(true) :
              intersection(false);
            if (temp.ray_intersects(curr_point,
                                    prev_point,
                                    ray,
                                    target_density,
                                    epsilon)) {
              temp.geo_div_id = gd_id;
              intersections.push_back(temp);
            }
            prev_point.x = curr_point.x;
            prev_point.y = curr_point.y;
          }
}

std::vector<std::vector<intersection> >
 InsetState::scanlines_parallel_to_axis(char grid_side, unsigned int res) const
{
  unsigned int grid_length = (grid_side == 'x') ? ly() : lx();
  unsigned int n_rays = grid_length * res;
  std::vector<std::vector<intersection> > scanlines(n_rays);

  // Iterate through GeoDivs in inset_state
  for (auto gd : geo_divs()) {

    // Find target density
    double target_density = target_areas_.at(gd.id()) / gd.area();

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
        // y = k and y = k+1 or x = k and x = k + 1, depending on grid_side
        for (double ray = k + (1.0/res)/2;
             ray < k + 1;
             ray += (1.0/res)) {

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
          // without entering or exiting the polygon. The case is the same when
          // x = ray_x goes exactly through curr_point.
          double epsilon = 1e-6 / res;

          // First we run the algorithm on the exterior ring
          add_intersections(intersections,
                            pwh.outer_boundary(),
                            ray,
                            target_density,
                            epsilon,
                            gd.id(),
                            grid_side);

          // Run algorithm on each hole
          for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
            add_intersections(intersections,
                              (*h),
                              ray,
                              target_density,
                              epsilon,
                              gd.id(),
                              grid_side);
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

// Returns line segments highlighting intersection points using scans above
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
