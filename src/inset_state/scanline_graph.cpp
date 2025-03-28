#include "inset_state.hpp"

// TODO: THE OUTPUT FROM intersec_with_parallel_to() ALWAYS
// SEEM TO COME WITH A NEED TO SORT AFTERWARDS. SHOULD SORTING BECOME PART of
// intersec_with_parallel_to() TO SAVE TYPING ELSEWHERE?

std::vector<std::vector<intersection> > InsetState::intersec_with_parallel_to(
  char axis,
  unsigned int resolution) const
{
  if (axis != 'x' && axis != 'y') {
    std::cerr << "ERROR: Invalid axis in " << __func__ << "()" << std::endl;
    exit(984320);
  }
  const unsigned int grid_length = (axis == 'x' ? ly_ : lx_);
  const unsigned int n_rays = grid_length * resolution;
  std::vector<std::vector<intersection> > scanlines(n_rays);

  // Iterate over GeoDivs in inset_state
  for (const auto &gd : geo_divs_) {

    // Find target density
    const double target_density = target_area_at(gd.id()) / gd.area();

    // Iterate over "polygons with holes" in inset_state
    for (size_t i = 0; const auto &pwh : gd.polygons_with_holes()) {
      const Bbox bb = pwh.bbox();
      double min_lim, max_lim;
      if (axis == 'x') {
        min_lim = bb.ymin();
        max_lim = bb.ymax();
      } else {
        min_lim = bb.xmin();
        max_lim = bb.xmax();
      }
      // Iterate over coordinates in bounding box of pwh
      // TODO: In the past, switching to the commented lines created problems
      // This should not be the case, as if floor and ceil return too low or
      // too high values, the rays will not intersect the polygon.
      // Investigate why the commented lines caused problems in the past.
      for (int k = min_lim - 1;
      // for (int k = floor(min_lim) - 1;
      // for (int k = min_lim;
      // for (int k = floor(min_lim);
      // for (unsigned int k = static_cast<unsigned int>(std::max(0, (int)floor(min_lim) - 1));
           k <= ceil(max_lim) + 1;
           // Same issue as above TODO here.
          //  k <= std::min((unsigned int)(max_lim), grid_length - 1);
           ++k) {

        // If the rays are in x-direction, iterate over each ray between the
        // grid lines y = k and y = k+1. Otherwise, iterate over each
        // ray between x = k and x = k+1.
        for (double ray = k + 0.5 / resolution; ray < k + 1;
             ray += (1.0 / resolution)) {

          // Temporary vector of intersections for this particular ray
          std::vector<intersection> intersections;

          // The following algorithm works by iterating over "resolution" rays
          // in each cell. // For each ray, we iterate over each edge in a
          // polygon and store all intersections. Afterward, we add the
          // appropriate densities. We add a small value `epsilon` to the ray
          // coordinate so that we assign correct densities if the ray with
          // equation y = ray_y or x = ray_x goes exactly through curr_point.
          // The addition ensures that, if there is any intersection, it is
          // only counted once. This method also correctly detects whether the
          // ray touches the point without entering or exiting the polygon.
          const double epsilon = 1e-6 / resolution;

          // Run algorithm on exterior ring
          add_intersections(
            intersections,
            pwh.outer_boundary(),
            ray,
            target_density,
            epsilon,
            gd.id(),
            i,
            axis);

          // Run algorithm on each hole
          for (const auto &h : pwh.holes()) {
            add_intersections(
              intersections,
              h,
              ray,
              target_density,
              epsilon,
              gd.id(),
              i,
              axis);
          }

          // Check whether the number of intersections is odd
          if (intersections.size() % 2 != 0) {
            std::cerr << "ERROR: Incorrect Topology. "
                      << "Number of intersections: " << intersections.size()
                      << ". "
                      << axis << "-coordinate: " << ray << ". "
                      << "Intersection points: ";

            for (auto &intersection : intersections) {
              std::cerr << (axis == 'x' ? intersection.x() : intersection.y())
                        << " ";
            }
            std::cerr << std::endl << std::endl;
            _Exit(932875);
          }
          std::sort(intersections.begin(), intersections.end());

          // Assign directions to intersections and add sorted vector of
          // intersections to vector `scanlines`
          const auto index = static_cast<unsigned int>(
            round((ray - 0.5 / resolution) * resolution));
          for (unsigned int l = 0; l < intersections.size(); ++l) {
            intersections[l].ray_enters = (l % 2 == 0);
            scanlines[index].push_back(intersections[l]);
          }
        }
      }
    }
  }
  return scanlines;
}

// Returns line segments highlighting intersection points using scans above
std::vector<Segment> InsetState::intersecting_segments(
  unsigned int resolution) const
{
  std::vector<Segment> int_segments;
  for (char axis : {'x', 'y'}) {
    const std::vector<std::vector<intersection> > scanlines =
      intersec_with_parallel_to(axis, resolution);
    const unsigned int grid_length = (axis == 'x' ? ly_ : lx_);

    // Iterate over rows (if axis is 'x') or columns
    for (unsigned int k = 0; k < grid_length; ++k) {

      // Iterate over rays in this row or column
      for (double ray = k + 0.5 / resolution; ray < k + 1;
           ray += (1.0 / resolution)) {

        // Intersections for one ray
        std::vector<intersection> intersec =
          scanlines[static_cast<unsigned int>(
            round((ray - (0.5 / resolution) * resolution)))];

        // Sort intersections in ascending order
        std::sort(intersec.begin(), intersec.end());

        // Check whether intersection enters twice or exits twice
        for (size_t l = 0; l + 1 < intersec.size(); ++l) {
          if (
            intersec[l].ray_enters == intersec[l + 1].ray_enters &&
            intersec[l].ray_enters && l + 2 <= intersec.size()) {
            Segment temp;
            if (axis == 'x' && intersec[l + 1].x() != intersec[l + 2].x()) {
              temp = Segment(
                Point(intersec[l + 1].x(), ray),
                Point(intersec[l + 2].x(), ray));
              int_segments.push_back(temp);
            } else if (intersec[l + 1].y() != intersec[l + 2].y()) {
              temp = Segment(
                Point(ray, intersec[l + 1].y()),
                Point(ray, intersec[l + 2].y()));
              int_segments.push_back(temp);
            }
          }
        }
      }
    }
  }
  return int_segments;
}
