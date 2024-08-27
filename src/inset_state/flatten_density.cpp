#include "constants.hpp"
#include "inset_state.hpp"
#include "interpolate_bilinearly.hpp"

// Function to calculate the velocity at the grid points (x, y) with x =
// 0.5, 1.5, ..., lx-0.5 and y = 0.5, 1.5, ..., ly-0.5 at time t
void calculate_velocity(
  double t,
  const FTReal2d &grid_fluxx_init,
  const FTReal2d &grid_fluxy_init,
  const FTReal2d &rho_ft,
  const FTReal2d &rho_init,
  boost::multi_array<double, 2> &grid_vx,
  boost::multi_array<double, 2> &grid_vy,
  const unsigned int lx,
  const unsigned int ly)
{
#pragma omp parallel for default(none) shared( \
    grid_fluxx_init,                           \
      grid_fluxy_init,                         \
      grid_vx,                                 \
      grid_vy,                                 \
      lx,                                      \
      ly,                                      \
      rho_ft,                                  \
      rho_init,                                \
      t)
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      double rho = rho_ft(0, 0) + (1.0 - t) * (rho_init(i, j) - rho_ft(0, 0));
      grid_vx[i][j] = -grid_fluxx_init(i, j) / rho;
      grid_vy[i][j] = -grid_fluxy_init(i, j) / rho;
    }
  }
}

bool all_points_are_in_domain(
  double delta_t,
  const boost::multi_array<Point, 2> &proj,
  const boost::multi_array<Vector, 2> &v_intp,
  const unsigned int lx,
  const unsigned int ly)
{
  // Return false if and only if there exists a point that would be outside
  // [0, lx] x [0, ly]
  bool in_domain = true;

#pragma omp parallel for reduction(&& : in_domain) default(none) \
  shared(delta_t, proj, v_intp, lx, ly)
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      double x = (proj)[i][j].x() + 0.5 * delta_t * (v_intp)[i][j].x();
      double y = (proj)[i][j].y() + 0.5 * delta_t * (v_intp)[i][j].y();
      if (x < 0.0 || x > lx || y < 0.0 || y > ly) {
        in_domain = false;
      }
    }
  }
  return in_domain;
}

// Function to integrate the equations of motion with the fast flow-based
// method
void InsetState::flatten_density()
{
  // Constants for the numerical integrator
  const double inc_after_acc = 1.1;
  const double dec_after_not_acc = 0.75;
  const double abs_tol = (std::min(lx_, ly_) * 1e-6);

  // Resize proj_ multi-array if running for the first time
  if (proj_.shape()[0] != lx_ || proj_.shape()[1] != ly_) {
    proj_.resize(boost::extents[lx_][ly_]);
  }

#pragma omp parallel for default(none) shared(proj_)
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      proj_[i][j] = Point(i + 0.5, j + 0.5);
    }
  }

  // Allocate memory for the velocity grid
  boost::multi_array<double, 2> grid_vx(boost::extents[lx_][ly_]);
  boost::multi_array<double, 2> grid_vy(boost::extents[lx_][ly_]);

  // eul[i][j] will be the new position of proj_[i][j] proposed by a simple
  // Euler step: move a full time interval delta_t with the velocity at time t
  // and position (proj_[i][j].x, proj_[i][j].y)
  boost::multi_array<Point, 2> eul(boost::extents[lx_][ly_]);

  // mid[i][j] will be the new displacement proposed by the midpoint
  // method (see comment below for the formula)
  boost::multi_array<Point, 2> mid(boost::extents[lx_][ly_]);

  // (vx_intp, vy_intp) will be the velocity at position (proj_.x, proj_.y) at
  // time t
  boost::multi_array<Vector, 2> v_intp(boost::extents[lx_][ly_]);

  // (vx_intp_half, vy_intp_half) will be the velocity at the midpoint
  // (proj_.x + 0.5*delta_t*vx_intp, proj_.y + 0.5*delta_t*vy_intp) at time
  // t + 0.5*delta_t
  boost::multi_array<Vector, 2> v_intp_half(boost::extents[lx_][ly_]);

  // Initialize the Fourier transforms of gridvx[] and gridvy[] at
  // every point on the lx_-times-ly_ grid at t = 0. We must typecast lx_ and
  // ly_ as double-precision numbers. Otherwise, the ratios in the denominator
  // will evaluate as zero.
  const double dlx = lx_;
  const double dly = ly_;

  // We temporarily insert the Fourier coefficients for the x-components and
  // y-components of the flux vector into grid_fluxx_init and grid_fluxy_init.
  // The reason for `+1` in `di+1` stems from the RODFT10 formula at:
  // https://www.fftw.org/fftw3_doc/1d-Real_002dodd-DFTs-_0028DSTs_0029.html
#pragma omp parallel for default(none) shared(dlx, dly)
  for (unsigned int i = 0; i < lx_ - 1; ++i) {
    double di = i;
    for (unsigned int j = 0; j < ly_; ++j) {
      double denom =
        pi * ((di + 1) / dlx + (j / (di + 1)) * (j / dly) * (dlx / dly));
      grid_fluxx_init_(i, j) = -rho_ft_(i + 1, j) / denom;
    }
  }
  for (unsigned int j = 0; j < ly_; ++j) {
    grid_fluxx_init_(lx_ - 1, j) = 0.0;
  }
#pragma omp parallel for default(none) shared(dlx, dly)
  for (unsigned int i = 0; i < lx_; ++i) {
    double di = i;
    for (unsigned int j = 0; j < ly_ - 1; ++j) {
      double denom =
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
  double delta_t = 1e-2;  // Initial time step.
  unsigned int iter = 0;
  const unsigned int max_iter = 300;

  // Integrate
  while (t < 1.0 && iter <= max_iter) {
    calculate_velocity(
      t,
      grid_fluxx_init_,
      grid_fluxy_init_,
      rho_ft_,
      rho_init_,
      grid_vx,
      grid_vy,
      lx_,
      ly_);
#pragma omp parallel for default(none) shared(proj_, v_intp, grid_vx, grid_vy)
    for (unsigned int i = 0; i < lx_; ++i) {
      for (unsigned int j = 0; j < ly_; ++j) {
        // We know, either because of the initialization or because of the
        // check at the end of the last iteration, that (proj_.x, proj_.y)
        // is inside the rectangle [0, lx_] x [0, ly_]. This fact guarantees
        // that interpolate_bilinearly() is given a point that cannot cause it
        // to fail.
        double vx = interpolate_bilinearly(
          proj_[i][j].x(),
          proj_[i][j].y(),
          grid_vx,
          'x',
          lx_,
          ly_);
        double vy = interpolate_bilinearly(
          proj_[i][j].x(),
          proj_[i][j].y(),
          grid_vy,
          'y',
          lx_,
          ly_);
        v_intp[i][j] = Vector(vx, vy);
      }
    }

    bool accept = false;
    while (!accept) {
// Simple Euler step
#pragma omp parallel for default(none) shared(proj_, v_intp, delta_t, eul)
      for (unsigned int i = 0; i < lx_; ++i) {
        for (unsigned int j = 0; j < ly_; ++j) {
          eul[i][j] = Point(
            proj_[i][j].x() + v_intp[i][j].x() * delta_t,
            proj_[i][j].y() + v_intp[i][j].y() * delta_t);
        }
      }

      // Use "explicit midpoint method"
      // x <- x + delta_t * v_x(x + 0.5*delta_t*v_x(x,y,t),
      //                        y + 0.5*delta_t*v_y(x,y,t),
      //                        t + 0.5*delta_t)
      // and similarly for y.
      calculate_velocity(
        t + 0.5 * delta_t,
        grid_fluxx_init_,
        grid_fluxy_init_,
        rho_ft_,
        rho_init_,
        grid_vx,
        grid_vy,
        lx_,
        ly_);

      // Make sure we do not pass a point outside [0, lx_] x [0, ly_] to
      // interpolate_bilinearly(). Otherwise, decrease the time step below and
      // try again.
      accept = all_points_are_in_domain(delta_t, proj_, v_intp, lx_, ly_);
      if (accept) {
        // Okay, we can run interpolate_bilinearly()
#pragma omp parallel for default(none) shared( \
    abs_tol,                                   \
      accept,                                  \
      delta_t,                                 \
      eul,                                     \
      grid_vx,                                 \
      grid_vy,                                 \
      mid,                                     \
      v_intp,                                  \
      v_intp_half)
        for (unsigned int i = 0; i < lx_; ++i) {
          for (unsigned int j = 0; j < ly_; ++j) {
            double vx_half = interpolate_bilinearly(
              proj_[i][j].x() + 0.5 * delta_t * v_intp[i][j].x(),
              proj_[i][j].y() + 0.5 * delta_t * v_intp[i][j].y(),
              grid_vx,
              'x',
              lx_,
              ly_);
            double vy_half = interpolate_bilinearly(
              proj_[i][j].x() + 0.5 * delta_t * v_intp[i][j].x(),
              proj_[i][j].y() + 0.5 * delta_t * v_intp[i][j].y(),
              grid_vy,
              'y',
              lx_,
              ly_);
            v_intp_half[i][j] = Vector(vx_half, vy_half);
            mid[i][j] = Point(
              proj_[i][j].x() + v_intp_half[i][j].x() * delta_t,
              proj_[i][j].y() + v_intp_half[i][j].y() * delta_t);

            // Do not accept the integration step if the maximum squared
            // difference between the Euler and midpoint proposals exceeds
            // abs_tol. Neither should we accept the integration step if one
            // of the positions wandered out of the domain. If one of these
            // problems occurred, decrease the time step.
            double sq_dist = CGAL::squared_distance(mid[i][j], eul[i][j]);
            if (
              sq_dist > abs_tol || mid[i][j].x() < 0.0 ||
              mid[i][j].x() > lx_ || mid[i][j].y() < 0.0 ||
              mid[i][j].y() > ly_) {
              accept = false;
            }
          }
        }
      }
      if (!accept) {
        delta_t *= dec_after_not_acc;
      }
    }

    // Control output and update for next iteration
    if (iter % 10 == 0) {
      std::cerr << "iter = " << iter << ", t = " << t
                << ", delta_t = " << delta_t << "\n";
    }

    // When we get here, the integration step was accepted
    t += delta_t;
    ++iter;
    proj_ = mid;
    delta_t *= inc_after_acc;  // Try a larger step next time
  }
}

bool all_map_points_are_in_domain(
  const double delta_t,
  const std::unordered_map<Point, Point> &proj_map,
  const std::unordered_map<Point, Vector> &v_intp,
  const unsigned int lx,
  const unsigned int ly)
{
  // Return false if and only if there exists a point that would be outside
  // [0, lx] x [0, ly]
  for (const auto &[key, val] : proj_map) {
    double x = val.x() + 0.5 * delta_t * v_intp.at(key).x();
    double y = val.y() + 0.5 * delta_t * v_intp.at(key).y();

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

double calculate_velocity_for_point(
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

bool InsetState::flatten_density_with_node_vertices()
{
  std::cerr << "In flatten_density_with_node_vertices()" << std::endl;

  // Constants for the numerical integrator
  const double inc_after_acc = 1.5;
  const double dec_after_not_acc = 0.5;
  const double abs_tol = (std::min(lx_, ly_) * 1e-6);

  // Clear previous triangle transformation data
  proj_qd_.triangle_transformation.clear();

  for (const Point &pt : unique_quadtree_corners_) {
    proj_qd_.triangle_transformation.insert_or_assign(pt, pt);
  }

  // eul[(i, j)] will be the new position of
  // proj_qd_.triangle_transformation[(i, j)] proposed by a simple Euler step:
  // move a full time interval delta_t with the velocity at time t and position
  // (proj_qd_.triangle_transformation[(i, j)].x,
  // proj_qd_.triangle_transformation[(i, j)].y)
  std::unordered_map<Point, Point> eul;

  // mid[(i, j)] will be the new displacement proposed by the midpoint
  // method (see comment below for the formula)
  std::unordered_map<Point, Point> mid;

  // v_intp[(i, j)] will be the velocity at position
  // (proj_qd_.triangle_transformation[(i, j)].x,
  // proj_qd_.triangle_transformation[(i, j)].y) at time t
  std::unordered_map<Point, Vector> v_intp;

  // v_intp_half[(i, j)] will be the velocity at the midpoint
  // (proj_qd_.triangle_transformation[(i, j)].x + 0.5 * delta_t * v_intp[(i,
  // j)].x, proj_qd_.triangle_transformation[(i, j)].y + 0.5 * delta_t *
  // v_intp[(i, j)].y) at time t + 0.5 * delta_t
  std::unordered_map<Point, Vector> v_intp_half;

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
    std::function<double(unsigned int, unsigned int, char)>
      cal_velocity_at_current_time =
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

    for (const auto &[key, val] : proj_qd_.triangle_transformation) {

      // We know, either because of the initialization or because of the
      // check at the end of the last iteration, that
      // proj_qd_.triangle_transformation[(i, j)] is inside the rectangle
      // [0, lx_] x [0, ly_]. This fact guarantees that
      // interpolate_bilinearly() is given a point that cannot cause it to
      // fail.
      Vector v_intp_val(
        interpolate_bilinearly(
          val.x(),
          val.y(),
          cal_velocity_at_current_time,
          'x',
          lx_,
          ly_),
        interpolate_bilinearly(
          val.x(),
          val.y(),
          cal_velocity_at_current_time,
          'y',
          lx_,
          ly_));
      v_intp.insert_or_assign(key, v_intp_val);
    }

    bool accept = false;
    while (!accept) {

      // Simple Euler step.
      for (const auto &[key, val] : proj_qd_.triangle_transformation) {
        Point eul_val(
          val.x() + v_intp[key].x() * delta_t,
          val.y() + v_intp[key].y() * delta_t);
        eul.insert_or_assign(key, eul_val);
      }

      // Use "explicit midpoint method"
      // x <- x + delta_t * v_x(x + 0.5 * delta_t * v_x(x, y, t),
      //                        y + 0.5 * delta_t *v_y(x, y, t),
      //                        t + 0.5 * delta_t)
      // and similarly for y.
      std::function<double(unsigned int, unsigned int, char)>
        cal_velocity_at_mid_time =
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
      accept = all_map_points_are_in_domain(
        delta_t,
        proj_qd_.triangle_transformation,
        v_intp,
        lx_,
        ly_);

      // Check whether any Delaunay triangle projections flip the projected
      // triangle. This is an issue as it can lead to intersection during
      // project. We resolve by increasing the blur width and running the steps
      // again.
      for (Delaunay::Finite_faces_iterator fit =
             proj_qd_.dt.finite_faces_begin();
           fit != proj_qd_.dt.finite_faces_end();
           ++fit) {
        const Point p1_ori = fit->vertex(0)->point();
        const Point p2_ori = fit->vertex(1)->point();
        const Point p3_ori = fit->vertex(2)->point();

        const Point p1_proj = proj_qd_.triangle_transformation[p1_ori];
        const Point p2_proj = proj_qd_.triangle_transformation[p2_ori];
        const Point p3_proj = proj_qd_.triangle_transformation[p3_ori];

        const double area_ori_triangle = CGAL::area(p1_ori, p2_ori, p3_ori);
        const double area_proj_triangle =
          CGAL::area(p1_proj, p2_proj, p3_proj);

        // If the area of the original triangle and the area of the projected
        // triangle have different signs, then the triangle has flipped.
        if (area_ori_triangle * area_proj_triangle < 0) {
          std::cerr << "Delaunay triangle flipped detected. Increasing blur "
                       "width and running again."
                    << std::endl;
          return 0;
        }
      }
      if (accept) {

        // Okay, we can run interpolate_bilinearly()
        for (const auto &[key, val] : proj_qd_.triangle_transformation) {
          Vector v_intp_half_val(
            interpolate_bilinearly(
              val.x() + 0.5 * delta_t * v_intp[key].x(),
              val.y() + 0.5 * delta_t * v_intp[key].y(),
              cal_velocity_at_mid_time,
              'x',
              lx_,
              ly_),
            interpolate_bilinearly(
              val.x() + 0.5 * delta_t * v_intp[key].x(),
              val.y() + 0.5 * delta_t * v_intp[key].y(),
              cal_velocity_at_mid_time,
              'y',
              lx_,
              ly_));
          v_intp_half.insert_or_assign(key, v_intp_half_val);
          Point mid_val(
            val.x() + v_intp_half[key].x() * delta_t,
            val.y() + v_intp_half[key].y() * delta_t);
          mid.insert_or_assign(key, mid_val);

          // Do not accept the integration step if the maximum squared
          // difference between the Euler and midpoint proposals exceeds
          // abs_tol. Neither should we accept the integration step if one
          // of the positions wandered out of the domain. If one of these
          // problems occurred, decrease the time step.
          const double sq_dist = CGAL::squared_distance(mid[key], eul[key]);
          if (
            sq_dist > abs_tol || mid[key].x() < 0.0 || mid[key].x() > lx_ ||
            mid[key].y() < 0.0 || mid[key].y() > ly_) {
            accept = false;
          }
        }
      }
      if (!accept) {
        delta_t *= dec_after_not_acc;
      }
    }

    // Control output
    if (iter % 10 == 0) {
      std::cerr << "iter = " << iter << ", t = " << t
                << ", delta_t = " << delta_t << "\n";
    }

    // When we get here, the integration step was accepted
    t += delta_t;
    ++iter;

    // Update the triangle transformation map
    proj_qd_.triangle_transformation = mid;
    delta_t *= inc_after_acc;  // Try a larger step next time
  }

  // Return 1 if the integration was successful
  return 1;
}