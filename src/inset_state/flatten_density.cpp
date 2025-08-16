#include "constants.hpp"
#include "inset_state.hpp"
#include "interpolate_bilinearly.hpp"

bool InsetState::flatten_density()
{

  // Create Delaunay triangulation based on quadtree corners and plot
  create_and_refine_quadtree();

  if (args_.plot_quadtree) {
    write_quadtree(file_prefix_ + "_quadtree");
    write_delaunay_triangles(file_prefix_ + "a_delaunay_t", false);
  }

  if (!flatten_density_on_node_vertices()) {

    // Flatten density has failed. Increase blur width and try again
    timer.stop("Flatten Density");
    return false;
  }

  // Flatten density passed.
  return true;
}

static bool all_map_points_are_in_domain(
  const double delta_t,
  const std::vector<Point> &projection,
  const std::vector<Vector> &velocity,
  const unsigned int lx,
  const unsigned int ly)
{
  // Return false if and only if there exists a point that would be outside
  // [0, lx] x [0, ly]
  for (size_t i = 0; i < projection.size(); i++) {
    const auto &pos = projection[i];
    const auto &velo = velocity[i];
    double x = pos.x() + 0.5 * delta_t * velo.x();
    double y = pos.y() + 0.5 * delta_t * velo.y();

    // if close to 0 using EPS, make 0, or greater than lx or ly, make lx or ly
    if (abs(x) < dbl_epsilon || abs(x - lx) < dbl_epsilon)
      x = (abs(x) < dbl_epsilon) ? 0 : lx;
    if (abs(y) < dbl_epsilon || abs(y - ly) < dbl_epsilon)
      y = (abs(y) < dbl_epsilon) ? 0 : ly;

    if (x < 0.0 || x > lx || y < 0.0 || y > ly) {
      return false;
    }
  }
  return true;
}

static inline double calculate_velocity_for_point(
  unsigned int i,
  unsigned int j,
  char direction,
  double t,
  const FTReal2d &grid_fluxx_init,
  const FTReal2d &grid_fluxy_init,
  const FTReal2d &rho_ft,
  const FTReal2d &rho_init)
{
  double rho = rho_ft(0, 0) + (1.0 - t) * (rho_init(i, j) - rho_ft(0, 0));
  return (direction == 'x') ? (-grid_fluxx_init(i, j) / rho)
                            : (-grid_fluxy_init(i, j) / rho);
}

bool InsetState::flatten_density_on_node_vertices()
{
  timer.start("Flatten Density");
  std::cerr << "In flatten_density_on_node_vertices()" << std::endl;

  // Constants for the numerical integrator
  const double inc_after_acc = 1.5;
  const double dec_after_not_acc = 0.5;
  const double reject_delta_t_threshold = 1e-4;
  const double abs_tol = (std::min(lx_, ly_) * 1e-6);

  const size_t num_quadtree_corners = unique_quadtree_corners_.size();

  std::vector<Point> &projection = proj_data_.get_projection();

  projection.resize(num_quadtree_corners);

  for (size_t i = 0; i < num_quadtree_corners; ++i) {
    projection[i] = unique_quadtree_corners_[i];
  }

  // eul[(i, j)] will be the new position of
  // proj_qd_.triangle_transformation[(i, j)] proposed by a simple Euler step:
  // move a full time interval delta_t with the velocity at time t and position
  // (proj_qd_.triangle_transformation[(i, j)].x,
  // proj_qd_.triangle_transformation[(i, j)].y)
  std::vector<Point> eul(num_quadtree_corners);

  // mid[(i, j)] will be the new displacement proposed by the midpoint
  // method (see comment below for the formula)
  std::vector<Point> mid(num_quadtree_corners);

  // v_intp[(i, j)] will be the velocity at position
  // (proj_qd_.triangle_transformation[(i, j)].x,
  // proj_qd_.triangle_transformation[(i, j)].y) at time t
  std::vector<Vector> v_intp(num_quadtree_corners);

  // We must typecast lx_ and ly_ as double-precision numbers. Otherwise, the
  // ratios in the denominator will evaluate as zero.
  double dlx = lx_;
  double dly = ly_;

  // We temporarily insert the Fourier coefficients for the x-components and
  // y-components of the flux vector into grid_fluxx_init and grid_fluxy_init
  // The reason for +1 in di+1 stems from the RODFT10 formula at:
  // https://www.fftw.org/fftw3_doc/1d-Real_002dodd-DFTs-_0028DSTs_0029.html
#pragma omp parallel for default(none) shared(dlx, dly)
  for (unsigned int i = 0; i < lx_ - 1; ++i) {
    const double di = i;
    for (unsigned int j = 0; j < ly_; ++j) {
      const double denom =
        pi * ((di + 1) / dlx + (j / (di + 1)) * (j / dly) * (dlx / dly));
      grid_fluxx_init_(i, j) = -rho_ft_(i + 1, j) / denom;
    }
  }
  for (unsigned int j = 0; j < ly_; ++j) {
    grid_fluxx_init_(lx_ - 1, j) = 0.0;
  }
#pragma omp parallel for default(none) shared(dlx, dly)
  for (unsigned int i = 0; i < lx_; ++i) {
    const double di = i;
    for (unsigned int j = 0; j < ly_ - 1; ++j) {
      const double denom =
        pi * ((di / (j + 1)) * (di / dlx) * (dly / dlx) + (j + 1) / dly);
      grid_fluxy_init_(i, j) = -rho_ft_(i, j + 1) / denom;
    }
  }
  for (unsigned int i = 0; i < lx_; ++i) {
    grid_fluxy_init_(i, ly_ - 1) = 0.0;
  }

  // Compute the flux vector and store the result in grid_fluxx_init and
  // grid_fluxy_init
  execute_fftw_plans_for_flux();

  double t = 0.0;
  double delta_t = 0.30;  // Initial time step.
  unsigned int iter = 0;
  unsigned int max_iter = 300;

  // Integrate
  while (t < 1.0 && iter <= max_iter) {

    // calculate_velocity lambda function
    auto cal_velocity_at_current_time =
      [&](unsigned int i, unsigned int j, char direction) {
        return calculate_velocity_for_point(
          i,
          j,
          direction,
          t,
          grid_fluxx_init_,
          grid_fluxy_init_,
          rho_ft_,
          rho_init_);
      };

    for (size_t i = 0; i < num_quadtree_corners; ++i) {
      const auto &pos = projection[i];

      // We know, either because of the initialization or because of the
      // check at the end of the last iteration, that
      // proj_qd_.triangle_transformation[(i, j)] is inside the rectangle
      // [0, lx_] x [0, ly_]. This fact guarantees that
      // interpolate_bilinearly() is given a point that cannot cause it to
      // fail.
      Vector v_intp_val(
        interpolate_bilinearly(
          pos.x(),
          pos.y(),
          cal_velocity_at_current_time,
          'x',
          lx_,
          ly_),
        interpolate_bilinearly(
          pos.x(),
          pos.y(),
          cal_velocity_at_current_time,
          'y',
          lx_,
          ly_));
      v_intp[i] = v_intp_val;
    }

    bool accept = false;
    while (!accept) {

      // Simple Euler step.
      for (size_t i = 0; i < num_quadtree_corners; ++i) {
        const auto &pos = projection[i];
        Point eul_val(
          pos.x() + v_intp[i].x() * delta_t,
          pos.y() + v_intp[i].y() * delta_t);
        eul[i] = eul_val;
      }

      // Use "explicit midpoint method"
      // x <- x + delta_t * v_x(x + 0.5 * delta_t * v_x(x, y, t),
      //                        y + 0.5 * delta_t *v_y(x, y, t),
      //                        t + 0.5 * delta_t)
      // and similarly for y.
      auto cal_velocity_at_mid_time =
        [&](unsigned int i, unsigned int j, char direction) {
          return calculate_velocity_for_point(
            i,
            j,
            direction,
            t + 0.5 * delta_t,
            grid_fluxx_init_,
            grid_fluxy_init_,
            rho_ft_,
            rho_init_);
        };

      // Make sure we do not pass a point outside [0, lx_] x [0, ly_] to
      // interpolate_bilinearly(). Otherwise, decrease the time step below and
      // try again.
      accept =
        all_map_points_are_in_domain(delta_t, projection, v_intp, lx_, ly_);

      if (accept) {

        // Okay, we can run interpolate_bilinearly()
        for (size_t i = 0; i < num_quadtree_corners; ++i) {
          const auto &pos = projection[i];
          const auto &velo = v_intp[i];

          // v_intp_half[(i, j)] will be the velocity at the midpoint
          // (proj_qd_.triangle_transformation[(i, j)].x + 0.5 * delta_t *
          // v_intp[(i, j)].x, proj_qd_.triangle_transformation[(i, j)].y + 0.5
          // * delta_t * v_intp[(i, j)].y) at time t + 0.5 * delta_t
          Vector v_intp_half_val(
            interpolate_bilinearly(
              pos.x() + 0.5 * delta_t * velo.x(),
              pos.y() + 0.5 * delta_t * velo.y(),
              cal_velocity_at_mid_time,
              'x',
              lx_,
              ly_),
            interpolate_bilinearly(
              pos.x() + 0.5 * delta_t * velo.x(),
              pos.y() + 0.5 * delta_t * velo.y(),
              cal_velocity_at_mid_time,
              'y',
              lx_,
              ly_));

          mid[i] = Point(
            pos.x() + v_intp_half_val.x() * delta_t,
            pos.y() + v_intp_half_val.y() * delta_t);

          // Do not accept the integration step if the maximum squared
          // difference between the Euler and midpoint proposals exceeds
          // abs_tol. Neither should we accept the integration step if one
          // of the positions wandered out of the domain. If one of these
          // problems occurred, decrease the time step.
          const double sq_dist = CGAL::squared_distance(mid[i], eul[i]);
          if (
            sq_dist > abs_tol || mid[i].x() < 0.0 || mid[i].x() > lx_ ||
            mid[i].y() < 0.0 || mid[i].y() > ly_) {
            accept = false;
          }
        }
      }
      if (!accept) {
        delta_t *= dec_after_not_acc;
        if (delta_t < reject_delta_t_threshold) {
          std::cerr << "Time step became too small. Increasing blur width and "
                       "running again."
                    << std::endl;
          return false;
        }
      }
    }

    // Control output
    if (iter % 10 == 0) {
      std::cerr << "iter " << iter << ": t = " << t
                << ", delta_t = " << delta_t << "\n";
    }

    // When we get here, the integration step was accepted
    t += delta_t;
    ++iter;

    // Update the triangle transformation map
    std::swap(projection, mid);

    delta_t *= inc_after_acc;  // Try a larger step next time
  }

  // Return true if the integration was successful
  timer.stop("Flatten Density");
  return true;
}
