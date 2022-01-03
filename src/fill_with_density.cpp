#include "cartogram_info.h"
#include "write_eps.h"
#include "inset_state.h"

void InsetState::fill_with_density(bool plot_density)
{
  // Calculate the total current area and total target area. We assume that
  // target areas that were zero or missing in the input have already been
  // replaced by CartogramInfo::replace_missing_and_zero_target_areas().
  double total_current_area = 0.0;
  for (auto gd : this->geo_divs()) {
    total_current_area += gd.area();
  }
  double total_target_area = 0.0;
  for (auto gd : this->geo_divs()) {
    total_target_area += this->target_areas_at(gd.id());
  }
  double mean_density = total_target_area / total_current_area;
  FTReal2d rho_init = this->rho_init_;

  // Initially assign 0 to all densities
  for (unsigned int i = 0; i < this->lx(); ++i) {
    for (unsigned int j = 0; j < this->ly(); ++j) {
      rho_init(i, j) = 0;
    }
  }

  // Resolution with which we sample polygons. "res" is the number of
  // horizontal "test rays" between each of the ly consecutive horizontal
  // graticule lines.
  unsigned int res = default_res;

  // A vector (map_intersections) to store vectors of intersections
  int n_rays = static_cast<int>(this->ly() * res);
  std::vector<std::vector<intersection> > map_intersections(n_rays);

  // Density numerator and denominator for each graticule cell
  // A density of a graticule cell can be calculated with (rho_num / rho_den).
  // We initially assign 0 to all elements because we assume that all
  // graticule cells are outside any GeoDiv. Any graticule cell where rho_den
  // is 0 will get the mean_density
  std::vector<std::vector<double> >
  rho_num(this->lx(), std::vector<double> (this->ly(), 0));
  std::vector<std::vector<double> >
  rho_den(this->lx(), std::vector<double> (this->ly(), 0));

  // See horizontal_scans in scanline_graph.cpp for more information.
  map_intersections = this->horizontal_scans(res);

  // Determine rho's numerator and denominator:
  // - rho_num is the sum of (weight * target_density) for each segment of a
  //   ray that is inside a GeoDiv
  // - rho_den is the sum of the weights of a ray that is inside a GeoDiv.
  // The weight of a segment of a ray that is inside a GeoDiv is equal to
  // (length of the segment inside the geo_div) * (area error of the geodiv).
  // We cycle through y-coordinates in inset_state.
  for (unsigned int k = 0; k < this->ly(); ++k) {

    // Cycle through each of the "res" number of rays in one cell
    for (double ray_y = k + 0.5/res;
         ray_y < k + 1;
         ray_y += 1.0/res) {

      // Intersections for one ray
      std::vector<intersection> intersections =
        map_intersections[static_cast<int>(round((ray_y - 0.5/res) * res))];

      // Sort vector in ascending order of intersection
      std::sort(intersections.begin(), intersections.end());

      // If the ray has intersections, we fill any empty spaces between
      // GeoDivs. Please note that we cannot write the loop condition as:
      // l < intersections.size() - 1
      // because intersection.size() is an unsigned integer. If
      // intersections.size() equals zero, then the right-hand side would
      // evaluate to a large positive number instead of -1. In this case,
      // we would erroneously enter the loop.
      for (unsigned int l = 1; l + 1 < intersections.size(); l += 2) {
        double left_x = intersections[l].coord;
        double right_x = intersections[l + 1].coord;
        if (left_x != right_x) {
          if (ceil(left_x) == ceil(right_x)) {

            // The intersections are in the same graticule cell. The ray
            // enters and leaves a GeoDiv in this cell. We weigh the density
            // of the cell by the GeoDiv's area error.
            double weight =
              this->area_errors_at(intersections[l].geo_div_id) *
              (right_x - left_x);
            double target_dens = intersections[l].target_density;
            rho_num[ceil(left_x) - 1][k] += weight * target_dens;
            rho_den[ceil(left_x) - 1][k] += weight;
          }
        }

        // Fill last exiting intersection with GeoDiv where part of ray inside
        // the graticule cell is inside the GeoDiv
        unsigned int last_x = intersections.back().coord;
        double last_weight =
          this->area_errors_at(intersections.back().geo_div_id) *
          (ceil(last_x) - last_x);
        double last_target_density = intersections.back().target_density;
        rho_num[ceil(last_x) - 1][k] += last_weight * last_target_density;
        rho_den[ceil(last_x) - 1][k] += last_weight;
      }

      // Fill GeoDivs by iterating through intersections
      for (unsigned int l = 0; l < intersections.size(); l += 2) {
        double left_x = intersections[l].coord;
        double right_x = intersections[l + 1].coord;

        // Check for intersection of polygons, holes and GeoDivs
        // TODO: Decide whether to comment out? (probably not)
        if (intersections[l].direction == intersections[l + 1].direction) {

          // Highlight where intersection is present
          std::cerr << "\nInvalid Geometry!" << std::endl;
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
            this->area_errors_at(intersections[l].geo_div_id);
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

  // Fill rho_init with the ratio of rho_num to rho_den
  for (unsigned int i = 0; i < this->lx(); ++i) {
    for (unsigned int j = 0; j < this->ly(); ++j) {
      if (rho_den[i][j] == 0) {
        rho_init(i, j) = mean_density;
      } else {
        rho_init(i, j) = rho_num[i][j] / rho_den[i][j];
      }
    }
  }
  if (plot_density) {
    std::string file_name =
      this->inset_name() +
      "_unblurred_density_" +
      std::to_string(this->n_finished_integrations()) +
      ".eps";
    std::cerr << "Writing " << file_name << std::endl;
    write_density_to_eps(file_name, rho_init.as_1d_array(), this);
  }
  this->execute_fftw_fwd_plan();
  return;
}
