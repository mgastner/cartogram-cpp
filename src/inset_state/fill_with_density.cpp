#include "../cartogram_info.h"
#include "../write_eps.h"
#include "../inset_state.h"

void InsetState::fill_with_density(bool plot_density)
{
  // Calculate the total current area and total target area. We assume that
  // target areas that were zero or missing in the input have already been
  // replaced by CartogramInfo::replace_missing_and_zero_target_areas().
  double total_current_area = 0.0;
  for (const auto &gd : geo_divs_) {
    total_current_area += gd.area();
  }
  double total_target_area = 0.0;
  for (const auto &gd : geo_divs_) {
    total_target_area += target_areas_at(gd.id());
  }
  double mean_density = total_target_area / total_current_area;
  FTReal2d rho_init = rho_init_;

  // Initially assign 0 to all densities
  for (unsigned int i = 0; i < lx(); ++i) {
    for (unsigned int j = 0; j < ly(); ++j) {
      rho_init(i, j) = 0;
    }
  }

  // Resolution with which we sample polygons. "resolution" is the number of
  // horizontal "test rays" between each of the ly consecutive horizontal
  // graticule lines.
  unsigned int resolution = default_resolution;

  // A vector (map_intersections) to store vectors of intersections
  unsigned int n_rays = ly() * resolution;
  std::vector<std::vector<intersection> > map_intersections(n_rays);

  // Density numerator and denominator for each graticule cell. The density of
  // a graticule cell can be calculated with (rho_num / rho_den). We initially
  // assign 0 to all elements because we assume that all graticule cells are
  // outside any GeoDiv. Any graticule cell where rho_den is 0 will be filled
  // with the mean_density.
  std::vector<std::vector<double> >
    rho_num(lx(), std::vector<double> (ly(), 0));
  std::vector<std::vector<double> >
    rho_den(lx(), std::vector<double> (ly(), 0));

  // See scanlines_parallel_to_axis in scanline_graph.cpp for more information
  map_intersections = scanlines_parallel_to_axis('x', resolution);

  // Determine rho's numerator and denominator:
  // - rho_num is the sum of (weight * target_density) for each segment of a
  //   ray that is inside a GeoDiv
  // - rho_den is the sum of the weights of a ray that is inside a GeoDiv.
  // The weight of a segment of a ray that is inside a GeoDiv is equal to
  // (length of the segment inside the geo_div) * (area error of the geodiv).

  // Iterate over each of the "resolution" rays between the graticule lines
  // y = k and y = k+1
  for (unsigned int k = 0; k < ly(); ++k) {
    for (double ray_y = k + 0.5/resolution;
         ray_y < k + 1;
         ray_y += 1.0/resolution) {

      // Intersections for one ray
      std::vector<intersection> intersections =
        map_intersections[static_cast<int>(round((ray_y - 0.5/resolution) * resolution))];

      // Sort vector in ascending order of intersection
      std::sort(intersections.begin(), intersections.end());

      // If the ray has intersections, we fill any empty spaces between
      // GeoDivs. Please note that we cannot write the loop condition as:
      // ll < intersections.size() - 1
      // because intersection.size() is an unsigned integer. If
      // intersections.size() equals zero, then the right-hand side would
      // evaluate to a large positive number instead of -1. In this case,
      // we would erroneously enter the loop.
      for (unsigned int ll = 1; ll + 1 < intersections.size(); ll += 2) {
        double left_x = intersections[ll].x();
        double right_x = intersections[ll + 1].x();
        if (left_x != right_x) {
          if (ceil(left_x) == ceil(right_x)) {

            // The intersections are in the same graticule cell. The ray
            // enters and leaves a GeoDiv in this cell. We weigh the density
            // of the cell by the GeoDiv's area error.
            double weight =
              area_errors_at(intersections[ll].geo_div_id) *
              (right_x - left_x);
            double target_dens = intersections[ll].target_density;
            rho_num[ceil(left_x) - 1][k] += weight * target_dens;
            rho_den[ceil(left_x) - 1][k] += weight;
          }
        }

        // Fill last exiting intersection with GeoDiv where part of ray inside
        // the graticule cell is inside the GeoDiv
        unsigned int last_x = intersections.back().x();
        double last_weight =
          area_errors_at(intersections.back().geo_div_id) *
          (ceil(last_x) - last_x);
        double last_target_density = intersections.back().target_density;
        rho_num[ceil(last_x) - 1][k] += last_weight * last_target_density;
        rho_den[ceil(last_x) - 1][k] += last_weight;
      }

      // Fill GeoDivs by iterating through intersections
      for (unsigned int ll = 0; ll < intersections.size(); ll += 2) {
        double left_x = intersections[ll].x();
        double right_x = intersections[ll + 1].x();

        // Check for intersection of polygons, holes and GeoDivs
        // TODO: Decide whether to comment out? (probably not)
        if (intersections[ll].ray_enters == intersections[ll + 1].ray_enters) {

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
            area_errors_at(intersections[ll].geo_div_id);
          double target_dens = intersections[ll].target_density;
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
  for (unsigned int i = 0; i < lx(); ++i) {
    for (unsigned int j = 0; j < ly(); ++j) {
      if (rho_den[i][j] == 0) {
        rho_init(i, j) = mean_density;
      } else {
        rho_init(i, j) = rho_num[i][j] / rho_den[i][j];
      }
    }
  }
  if (plot_density) {
    std::string file_name =
      inset_name() +
      "_unblurred_density_" +
      std::to_string(n_finished_integrations()) +
      ".eps";
    std::cerr << "Writing " << file_name << std::endl;
    write_density_to_eps(file_name, rho_init.as_1d_array(), this);
  }
  execute_fftw_fwd_plan();
  return;
}
