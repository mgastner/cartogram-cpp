#include "constants.hpp"
#include "inset_state.hpp"

void InsetState::fill_with_density_rays()
{
  std::cerr << "Filling density using ray-shooting method" << std::endl;

  timer.start("Fill with Density (Ray Shooting Method)");
  // We assume that target areas that were zero or missing in the input have
  // already been replaced by
  // CartogramInfo::replace_missing_and_zero_target_areas().
  // We also correct for a drift in the total inset area, in case it is
  // present, by adjusting exterior_density. The idea is to treat the exterior
  // area as if it were a polygon with target area
  // C * (lx * ly - initial_area_) and current area
  // C * (lx * ly - total_inset_area). The constant prefactor C is the same
  // in both areas and, thus, cancels out when taking the ratio. Note that
  // missing, zero, and near-zero areas must already be filled with surrogate
  // target areas in the desired proportion to the other polygons. That is,
  // we must call cart_info.replace_missing_and_zero_target_areas() before
  // calling this function.
  double exterior_density =
    (lx_ * ly_ - initial_area_) / (lx_ * ly_ - total_inset_area());
  // We weight GeoDivs according to area errors, including the exterior.
  // area error is defined as | ((area / target_area) - 1) |
  // lx * ly - initial_area is thus, the target area, and
  // lx * ly - total_inset_area is thus, the current area
  // Thus, area / target_area of the exterior is the inverse
  // of the exterior_density.
  double ext_weight = std::abs((1.0 / exterior_density) - 1);

#pragma omp parallel for default(none)

  // Initially assign zero to all densities
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      rho_init_(i, j) = 0;
    }
  }

  // Density numerator and denominator for each grid cell. The density of
  // a grid cell can be calculated with (rho_num / rho_den). We initially
  // assign zero to all elements because we assume that all grid cells
  // are outside any GeoDiv. Any grid cell where rho_den is zero will be
  // filled with the exterior_density.

  boost::multi_array<double, 2> rho_num(boost::extents[lx_][ly_]);
  boost::multi_array<double, 2> rho_den(boost::extents[lx_][ly_]);
  std::fill_n(rho_num.data(), rho_num.num_elements(), 0.0);
  std::fill_n(rho_den.data(), rho_den.num_elements(), 0.0);

  // Initialize an unordered set with all Polygons with Holes, so that we can
  // check at the end if all GeoDivs have been processed
  std::map<std::string, std::set<size_t>> unseen_pwhs;
  for (const auto &gd : geo_divs_) {
    for (size_t i = 0; i < gd.polygons_with_holes().size(); ++i) {
      unseen_pwhs[gd.id()].insert(i);
    }
  }

  // Resolution with which we sample polygons. "resolution" is the number of
  // horizontal "test rays" between each of the ly consecutive horizontal
  // grid lines.
  // Ensure that the number of rays per grid length is at least the default
  // value specified by default_resolution in include/constants.hpp.
  // Additionally, confirm that there are a minimum of long_grid_length *
  // resolution rays along the longer side of the lx*ly grid.
  unsigned int long_grid_length = std::max(lx_, ly_);
  const unsigned int resolution =
    (long_grid_length > default_long_grid_length)
      ? static_cast<unsigned int>(
          (default_resolution * default_long_grid_length) *
          (1.0 / long_grid_length))
      : default_resolution;

  auto intersections_with_rays = intersec_with_parallel_to('x', resolution);

  // Determine total_interior_density to get the average density of the
  // interior later
  double total_interior_density = 0;

  // Determine rho's numerator and denominator:
  // - rho_num is the sum of (weight * target_density) for each segment of a
  //   ray that is inside a GeoDiv.
  // - rho_den is the sum of the weights of a ray that is inside a GeoDiv.
  // The weight of a segment of a ray that is inside a GeoDiv is equal to
  // (length of the segment inside the geo_div) * (area error of the geodiv).
#pragma omp parallel for default(none)                         \
  shared(intersections_with_rays, rho_den, rho_num, std::cerr) \
  firstprivate(resolution, exterior_density, ext_weight)       \
  reduction(+ : total_interior_density)
  for (unsigned int k = 0; k < ly_; ++k) {

    // Iterate over each of the rays between the grid lines y = k and
    // y = k+1
    for (double y = k + 0.5 / resolution; y < k + 1; y += 1.0 / resolution) {

      // Intersections for one ray
      auto intersections_at_y =
        intersections_with_rays[static_cast<std::size_t>(
          std::lround(resolution * y - 0.5))];

      // Sort intersections in ascending order
      std::sort(intersections_at_y.begin(), intersections_at_y.end());

      // If the ray has intersections, we fill cells between
      // the two intersections with the target density of the current GeoDiv.

      // Please note that we cannot write the loop condition as:
      // i < intersections.size() - 1
      // because intersection.size() is an unsigned integer. If
      // intersections.size() equals zero, then the right-hand side would
      // evaluate to a large positive number instead of -1. In this case,
      // we would erroneously enter the loop.
      for (unsigned int i = 0; i + 1 < intersections_at_y.size(); i += 2) {

        // Fill the left-most cell
        const double left_x = intersections_at_y[i].x();
        const double right_x = intersections_at_y[i + 1].x();
        const double target_dens = intersections_at_y[i].target_density;
        double weight = area_error_at(intersections_at_y[i].geo_div_id);

        total_interior_density += (right_x - left_x) * target_dens;

        // Check for intersection of polygons, holes and GeoDivs
        // TODO: Decide whether to comment out? (probably not)
        if (
          intersections_at_y[i].ray_enters ==
          intersections_at_y[i + 1].ray_enters) {

          // Highlight where intersection is present
          std::cerr << "\n ERROR: Invalid Geometry!" << std::endl;
          std::cerr << "Intersection of Polygons/Holes/Geodivs" << std::endl;
          std::cerr << "Y-coordinate: " << y << std::endl;
          std::cerr << "Left X-coordinate: " << left_x << std::endl;
          std::cerr << "Right X-coordinate: " << right_x << std::endl;
          std::exit(8026519);
        }

        // Ray is entering a GeoDiv from empty space
        double prev_x =
          (i == 0) ?
                   // It is the first intersection
            floor(left_x)
                   :
                   // prev_x could be within the same grid cell
            std::max(floor(left_x), intersections_at_y[i - 1].x());

        double left_x_empty_offset = left_x - prev_x;

        // The entire segment of GeoDiv may be inside one cell
        double x_after_left_x = std::min(ceil(left_x), right_x);
        double gd_length_to_left = x_after_left_x - left_x;
        // weight = 1;

        // Fill leftmost cell with GeoDiv where part of ray inside the grid
        // rho_num[floor(left_x)][k] += left_x_empty_offset * exterior_density
        // + gd_length_to_left * target_dens;
        rho_num[static_cast<unsigned int>(floor(left_x))][k] +=
          left_x_empty_offset * exterior_density * ext_weight +
          gd_length_to_left * target_dens * weight;
        // rho_den[floor(left_x)][k] += left_x_empty_offset +
        // gd_length_to_left;
        rho_den[static_cast<unsigned int>(floor(left_x))][k] +=
          left_x_empty_offset * ext_weight + weight * gd_length_to_left;

        // Fill all other cells where GeoDiv is fully covering grid cell
        for (unsigned int m = static_cast<unsigned int>(std::ceil(left_x));
             m < static_cast<unsigned int>(std::floor(right_x));
             ++m) {
          rho_num[m][k] += weight * target_dens;
          rho_den[m][k] += weight;
        }

        // Ray is exiting a GeoDiv to empty space
        double next_x =
          (i + 2 == intersections_at_y.size())
            ?
            // It is the last intersection
            ceil(right_x)
            :
            // next_x could be within the same grid cell
            std::min(ceil(right_x), intersections_at_y[i + 2].x());

        double right_x_empty_offset =
          !almost_equal(next_x, right_x)
            ?
            // There is space between two GeoDivs, within the same cell,
            // which will be handled in the next iteration upon ray entering
            // next GeoDiv
            0
            :
            // The ray is exiting the GeoDiv and is "empty" until the end of
            // the grid cell
            next_x - right_x;

        double gd_length_to_right =
          floor(right_x) <= left_x
            ?
            // This entire segment of GeoDiv is inside one cell
            // The earlier part of the code has already accounted for the
            // density of this segment.
            0
            : right_x - floor(right_x);

        // rho_num[floor(right_x)][k] += right_x_empty_offset *
        // exterior_density + gd_length_to_right * target_dens;
        rho_num[static_cast<unsigned int>(floor(right_x))][k] +=
          right_x_empty_offset * exterior_density * ext_weight +
          gd_length_to_right * target_dens * weight;
        // rho_den[floor(right_x)][k] += right_x_empty_offset *
        // exterior_density + gd_length_to_right;
        rho_den[static_cast<unsigned int>(floor(right_x))][k] +=
          right_x_empty_offset * exterior_density * ext_weight +
          gd_length_to_right * weight;
      }
    }
  }

  // The total_interior_density has been added resolution number of times per
  // grid cell
  total_interior_density /= resolution;
  dens_mean_ = total_interior_density / total_inset_area();

  // If any Polygons remain unseen, process them individually
  for (const auto &[gd_id, pwh_set] : unseen_pwhs) {

    // Skip iteration in case no Polygons with Holes of this GeoDiv remain
    // unseen
    if (pwh_set.empty())
      continue;

    // Get reference to GeoDiv
    const auto &gd = geo_div_at_id(gd_id);

    // Calculate weight and target density for the polygon
    const double weight = area_error_at(gd_id);
    const double target_dens = target_area_at(gd_id) / gd.area();

    // Skip iteration in case GeoDiv has already converged
    if (weight < 0.01)
      continue;

    // Iterate over all unseen Polygons with Holes of this GeoDiv to find
    for (const auto &pwh_idx : pwh_set) {

      const auto &pwh = gd.polygons_with_holes()[pwh_idx];

      // Get bounding box of the polygon to get the grid cell that contains it
      Bbox pwh_bbox = pwh.bbox();
      unsigned int grid_i = static_cast<unsigned int>(
        std::floor((pwh_bbox.xmin() + pwh_bbox.xmax()) / 2.0));
      unsigned int grid_j = static_cast<unsigned int>(
        std::floor((pwh_bbox.ymin() + pwh_bbox.ymax()) / 2.0));

      // TODO: This may attribute a very high weight, which would've been
      // otherwise reduced by the length of the segment inside the GeoDiv.
      // Investigate the best settings for this
      rho_num[grid_i][grid_j] += weight * target_dens;
      rho_den[grid_i][grid_j] += weight;
    }
  }

  // Fill rho_init with the ratio of rho_num to rho_den
#pragma omp parallel for default(none) shared(rho_den, rho_num, rho_init_) \
  firstprivate(exterior_density)
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      if (almost_equal(rho_den[i][j], 0.0)) {
        rho_init_(i, j) = exterior_density;
      } else {
        rho_init_(i, j) = rho_num[i][j] / rho_den[i][j];
      }
    }
  }

  // Determine range of densities
  auto [min_iter, max_iter] = std::minmax_element(
    rho_init_.as_1d_array(),
    rho_init_.as_1d_array() + lx_ * ly_);

  dens_min_ = *min_iter;
  exterior_density_ = exterior_density;
  dens_max_ = *max_iter;

  execute_fftw_fwd_plan();
  timer.stop("Fill with Density (Ray Shooting Method)");
}
