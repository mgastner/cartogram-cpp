#include "../inset_state.h"
#include "../constants.h"

std::vector<std::vector<intersection> >
InsetState::intersections_with_rays_parallel_to_axis(
  char axis,
  unsigned int resolution) const
{
  if (axis != 'x' && axis != 'y') {
    std::cerr << "Invalid axis in intersections_with_rays_parallel_to_axis()"
              << std::endl;
    exit(984320);
  }
  unsigned int grid_length = (axis == 'x' ? ly_ : lx_);
  unsigned int n_rays = grid_length * resolution;
  std::vector<std::vector<intersection> > scanlines(n_rays);

  // Iterate over GeoDivs in inset_state
  for (const auto &gd : geo_divs_) {

    // Find target density
    double target_density = target_areas_.at(gd.id()) / gd.area();

    // Iterate over "polygons with holes" in inset_state
    for (const auto &pwh : gd.polygons_with_holes()) {
      Bbox bb = pwh.bbox();
      double min_lim, max_lim;
      if (axis == 'x') {
        min_lim = bb.ymin();
        max_lim = bb.ymax();
      } else {
        min_lim = bb.xmin();
        max_lim = bb.xmax();
      }

      // Iterate over coordinates in bounding box of pwh
      for (unsigned int k = floor(min_lim) - 1;
           k <= ceil(max_lim) + 1;
           ++k) {

        // If the rays are in x-direction, iterate over each ray between the
        // graticule lines y = k and y = k+1. Otherwise, iterate over each
        // ray between x = k and x = k+1.
        for (double ray = k + (1.0/resolution)/2;
             ray < k + 1;
             ray += (1.0/resolution)) {

          // Temporary vector of intersections for this particular ray
          std::vector<intersection> intersections;

          // The following algorithm works by iterating over "resolution" rays
          // in each cell. // For each ray, we iterate over each edge in a
          // polygon and store all intersections. Afterward, we add the
          // appropriate densities. We add a small value `epsilon` to the ray
          // coordinate so that we assign correct densities if the ray with
          // equation y = ray_y or x = ray_x goes exactly through curr_point.
          // The addition ensures that, if there is any intersection, it is only
          // counted once. This method also correctly detects whether the ray
          // touches the point without entering or exiting the polygon.
          double epsilon = 1e-6 / resolution;

          // Run algorithm on exterior ring
          add_intersections(intersections,
                            pwh.outer_boundary(),
                            ray,
                            target_density,
                            epsilon,
                            gd.id(),
                            axis);

          // Run algorithm on each hole
          for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
            add_intersections(intersections,
                              *h,
                              ray,
                              target_density,
                              epsilon,
                              gd.id(),
                              axis);
          }

          // Check whether the number of intersections is odd
          if (intersections.size() % 2 != 0) {
            std::cerr << "Incorrect Topology.\n"
                      << "Number of intersections: "
                      << intersections.size()
                      << "\n"
                      << axis
                      << "-coordinate: "
                      << ray
                      << "\n"
                      << "Intersection points: "
                      << std::endl;

            for (unsigned int l = 0; l < intersections.size(); ++l) {
              std::cerr
                << (axis == 'x' ? intersections[l].x() : intersections[l].y())
                << std::endl;
            }
            std::cerr << std::endl << std::endl;
            _Exit(932875);
          }
          std::sort(intersections.begin(), intersections.end());

          // Assign directions to intersections and add sorted vector of
          // intersections to vector `scanlines`
          int index = round((ray - 0.5/resolution) * resolution);
          for (unsigned int l = 0; l < intersections.size(); ++l) {
            intersections[l].ray_enters = (l%2 == 0);
            scanlines[index].push_back(intersections[l]);
          }
        }
      }
    }
  }
  return scanlines;
}

// Creates continuity/adjacency graph using horizontal and vertical scans above
void InsetState::create_contiguity_graph(unsigned int resolution)
{
  // Calculate horizontal and vertical scanlines
  for (char axis : {'x', 'y'}) {
    std::vector<std::vector<intersection> > scanlines =
      intersections_with_rays_parallel_to_axis(axis, resolution);
    unsigned int max_k = (axis == 'x' ? ly_ : lx_);

    // Iterate over rows (if is_x_axis) or columns
    for (unsigned int k = 0; k < max_k; ++k) {

      // Iterate over rays in this row or column
      for (double ray = k + (1.0/resolution)/2;
           ray < k + 1;
           ray += (1.0/resolution)) {

        // Intersections for one ray
        std::vector<intersection> intersections =
          scanlines[round(((ray - (1.0/resolution)/2.0) * resolution))];

        // Sort intersections in ascending order
        sort(intersections.begin(), intersections.end());
        int size = intersections.size() - 1;

        // Find adjacent GeoDivs by iterating over intersections
        for (int l = 1; l < size; l += 2) {
          double coord_1 = intersections[l].x();
          double coord_2 = intersections[l + 1].x();
          std::string gd_1 = intersections[l].geo_div_id;
          std::string gd_2 = intersections[l + 1].geo_div_id;

          // Update adjacency
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

// Returns line segments highlighting intersection points using scans above
const std::vector<Segment>
InsetState::intersecting_segments(unsigned int resolution) const
{
  std::vector<Segment> int_segments;
  for (char axis : {'x', 'y'}) {
    std::vector<std::vector<intersection> > scanlines =
      intersections_with_rays_parallel_to_axis(axis, resolution);
    unsigned int max_k = (axis == 'x' ? ly_ : lx_);

    // Iterate over rows (if is_x_axis) or columns
    for (unsigned int k = 0; k < max_k; ++k) {

      // Iterate over rays in this row or column
      for (double ray = k + (1.0/resolution)/2;
           ray < k + 1;
           ray += (1.0/resolution)) {

        // Intersections for one ray
        std::vector<intersection> intersections =
          scanlines[round(((ray - (1.0/resolution)/2.0) * resolution))];

        // Sort intersections in ascending order
        std::sort(intersections.begin(), intersections.end());

        // Check whether intersection enters twice or exits twice
        for (size_t l = 0; l + 1 < intersections.size(); ++l) {
          if (intersections[l].ray_enters == intersections[l + 1].ray_enters &&
              intersections[l].ray_enters &&
              l + 2 <= intersections.size()) {
            Segment temp;
            if (axis == 'x' &&
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
