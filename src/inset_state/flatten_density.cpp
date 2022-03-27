#include "constants.h"
#include "inset_state.h"
#include "interpolate_bilinearly.h"
#include <boost/multi_array.hpp>
#include <omp.h>

// Function to calculate the velocity at the grid points (x, y) with x =
// 0.5, 1.5, ..., lx-0.5 and y = 0.5, 1.5, ..., ly-0.5 at time t
void calculate_velocity(
    double t,
    FTReal2d &grid_fluxx_init,
    FTReal2d &grid_fluxy_init,
    FTReal2d &rho_ft,
    FTReal2d &rho_init_,
    boost::multi_array<double, 2> *grid_vx,
    boost::multi_array<double, 2> *grid_vy,
    const unsigned int lx,
    const unsigned int ly
) {
  double rho;

#pragma omp parallel for private(rho)
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      rho = rho_ft(0, 0) + (1.0 - t) * (rho_init_(i, j) - rho_ft(0,0));
      (*grid_vx)[i][j] = -grid_fluxx_init(i, j) / rho;
      (*grid_vy)[i][j] = -grid_fluxy_init(i, j) / rho;
    }
  }
  return;
}

bool all_points_are_in_domain(
    double delta_t,
    boost::multi_array<XYPoint, 2> *proj,
    boost::multi_array<XYPoint, 2> *v_intp,
    const unsigned int lx,
    const unsigned int ly
) {
  // Return false if and only if there exists a point that would be outside
  // [0, lx] x [0, ly]
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      double x = (*proj)[i][j].x + 0.5 * delta_t * (*v_intp)[i][j].x;
      double y = (*proj)[i][j].y + 0.5 * delta_t * (*v_intp)[i][j].y;
      if (x < 0.0 || x > lx || y < 0.0 || y > ly) {
        return false;
      }
    }
  }
  return true;
}

// Function to integrate the equations of motion with the fast flow-based
// method
void InsetState::flatten_density()
{
  std::cerr << "In flatten_density()" << std::endl;

  // Constants for the numerical integrator
  const double inc_after_acc = 1.1;
  const double dec_after_not_acc = 0.75;
  const double abs_tol = (std::min(lx_, ly_) * 1e-6);

  // Resize proj_ multi-array if running for the first time
  if (proj_.shape()[0] != lx_ || proj_.shape()[1] != ly_) {
    proj_.resize(boost::extents[lx_][ly_]);
  }
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      proj_[i][j].x = i + 0.5;
      proj_[i][j].y = j + 0.5;
    }
  }

  // Allocate memory for the velocity grid
  boost::multi_array<double, 2> grid_vx(boost::extents[lx_][ly_]);
  boost::multi_array<double, 2> grid_vy(boost::extents[lx_][ly_]);

  // Prepare Fourier transforms for the flux
  FTReal2d grid_fluxx_init;
  FTReal2d grid_fluxy_init;
  grid_fluxx_init.allocate(lx_, ly_);
  grid_fluxy_init.allocate(lx_, ly_);
  grid_fluxx_init.make_fftw_plan(FFTW_RODFT01, FFTW_REDFT01);
  grid_fluxy_init.make_fftw_plan(FFTW_REDFT01, FFTW_RODFT01);

  // eul[i][j] will be the new position of proj_[i][j] proposed by a simple
  // Euler step: move a full time interval delta_t with the velocity at time t
  // and position (proj_[i][j].x, proj_[i][j].y)
  boost::multi_array<XYPoint, 2> eul(boost::extents[lx_][ly_]);

  // mid[i][j] will be the new displacement proposed by the midpoint
  // method (see comment below for the formula)
  boost::multi_array<XYPoint, 2> mid(boost::extents[lx_][ly_]);

  // (vx_intp, vy_intp) will be the velocity at position (proj_.x, proj_.y) at
  // time t
  boost::multi_array<XYPoint, 2> v_intp(boost::extents[lx_][ly_]);

  // (vx_intp_half, vy_intp_half) will be the velocity at the midpoint
  // (proj_.x + 0.5*delta_t*vx_intp, proj_.y + 0.5*delta_t*vy_intp) at time
  // t + 0.5*delta_t
  boost::multi_array<XYPoint, 2> v_intp_half(boost::extents[lx_][ly_]);

  // Initialize the Fourier transforms of gridvx[] and gridvy[] at
  // every point on the lx_-times-ly_ grid at t = 0. We must typecast lx_ and ly_
  // as double-precision numbers. Otherwise the ratios in the denominator
  // will evaluate as zero.
  double dlx = lx_;
  double dly = ly_;

  // We temporarily insert the Fourier coefficients for the x-components and
  // y-components of the flux vector into grid_fluxx_init and grid_fluxy_init
  for (unsigned int i = 0; i < lx_-1; ++i) {
    double di = i;
    for (unsigned int j = 0; j < ly_; ++j) {
      double denom = pi * ((di+1)/dlx + (j/(di+1)) * (j/dly) * (dlx/dly));
      grid_fluxx_init(i, j) =
        -rho_ft_(i+1, j) / denom;
    }
  }
  for (unsigned int j = 0; j < ly_; ++j) {
    grid_fluxx_init(lx_-1, j) = 0.0;
  }
  for (unsigned int i=0; i<lx_; ++i) {
    double di = i;
    for (unsigned int j = 0; j < ly_-1; ++j) {
      double denom = pi * ((di/(j+1)) * (di/dlx) * (dly/dlx) + (j+1)/dly);
      grid_fluxy_init(i, j) =
        -rho_ft_(i, j+1) / denom;
    }
  }
  for (unsigned int i=0; i<lx_; ++i) {
    grid_fluxy_init(i, ly_-1) = 0.0;
  }

  // Compute the flux vector and store the result in grid_fluxx_init and
  // grid_fluxy_init
  grid_fluxx_init.execute_fftw_plan();
  grid_fluxy_init.execute_fftw_plan();
  double t = 0.0;
  double delta_t = 1e-2;  // Initial time step.
  unsigned int iter = 0;

  // Integrate
  while (t < 1.0) {
    calculate_velocity(
      t,
      grid_fluxx_init,
      grid_fluxy_init,
      rho_ft_,
      rho_init_,
      &grid_vx,
      &grid_vy,
      lx_,
      ly_
    );
#pragma omp parallel for
    for (unsigned int i = 0; i < lx_; ++i) {
      for (unsigned int j = 0; j < ly_; ++j) {

        // We know, either because of the initialization or because of the
        // check at the end of the last iteration, that (proj_.x, proj_.y)
        // is inside the rectangle [0, lx_] x [0, ly_]. This fact guarantees
        // that interpolate_bilinearly() is given a point that cannot cause it
        // to fail.
        v_intp[i][j].x = interpolate_bilinearly(
          proj_[i][j].x,
          proj_[i][j].y,
          &grid_vx,
          'x',
          lx_,
          ly_
        );
        v_intp[i][j].y = interpolate_bilinearly(
          proj_[i][j].x,
          proj_[i][j].y,
          &grid_vy,
          'y',
          lx_,
          ly_
        );
      }
    }
    bool accept = false;
    while (!accept) {

      // Simple Euler step.

#pragma omp parallel for
      for (unsigned int i = 0; i < lx_; ++i) {
        for (unsigned int j = 0; j < ly_; ++j) {
          eul[i][j].x = proj_[i][j].x + v_intp[i][j].x * delta_t;
          eul[i][j].y = proj_[i][j].y + v_intp[i][j].y * delta_t;
        }
      }

      // Use "explicit midpoint method"
      // x <- x + delta_t * v_x(x + 0.5*delta_t*v_x(x,y,t),
      //                        y + 0.5*delta_t*v_y(x,y,t),
      //                        t + 0.5*delta_t)
      // and similarly for y.
      calculate_velocity(
        t + 0.5*delta_t,
        grid_fluxx_init,
        grid_fluxy_init,
        rho_ft_,
        rho_init_,
        &grid_vx,
        &grid_vy,
        lx_,
        ly_
      );

      // Make sure we do not pass a point outside [0, lx_] x [0, ly_] to
      // interpolate_bilinearly(). Otherwise decrease the time step below and
      // try again.
      accept = all_points_are_in_domain(delta_t, &proj_, &v_intp, lx_, ly_);
      if (accept) {

        // Okay, we can run interpolate_bilinearly()

#pragma omp parallel for
        for (unsigned int i = 0; i < lx_; ++i) {
          for (unsigned int j = 0; j < ly_; ++j) {
            v_intp_half[i][j].x = interpolate_bilinearly(
              proj_[i][j].x + 0.5*delta_t*v_intp[i][j].x,
              proj_[i][j].y + 0.5*delta_t*v_intp[i][j].y,
              &grid_vx,
              'x',
              lx_,
              ly_
            );
            v_intp_half[i][j].y = interpolate_bilinearly(
              proj_[i][j].x + 0.5*delta_t*v_intp[i][j].x,
              proj_[i][j].y + 0.5*delta_t*v_intp[i][j].y,
              &grid_vy,
              'y',
              lx_,
              ly_
            );
            mid[i][j].x = proj_[i][j].x + v_intp_half[i][j].x * delta_t;
            mid[i][j].y = proj_[i][j].y + v_intp_half[i][j].y * delta_t;

            // Do not accept the integration step if the maximum squared
            // difference between the Euler and midpoint proposals exceeds
            // abs_tol. Neither should we accept the integration step if one
            // of the positions wandered out of the domain. If one of these
            // problems occurred, decrease the time step.
            const double sq_dist =
              (mid[i][j].x-eul[i][j].x) * (mid[i][j].x-eul[i][j].x)
              + (mid[i][j].y-eul[i][j].y) * (mid[i][j].y-eul[i][j].y);
            if (sq_dist > abs_tol ||
                mid[i][j].x < 0.0 || mid[i][j].x > lx_ ||
                mid[i][j].y < 0.0 || mid[i][j].y > ly_) {
              accept = false;
            }
          }
        }
      }
      if (!accept) {
        delta_t *= dec_after_not_acc;
      }
    }

    // Control ouput
    if (iter % 10 == 0) {
      std::cerr << "iter = "
                << iter
                << ", t = "
                << t
                << ", delta_t = "
                << delta_t
                << "\n";
    }

    // When we get here, the integration step was accepted
    t += delta_t;
    ++iter;
    proj_ = mid;
    delta_t *= inc_after_acc;  // Try a larger step next time
  }
  grid_fluxx_init.destroy_fftw_plan();
  grid_fluxy_init.destroy_fftw_plan();
  grid_fluxx_init.free();
  grid_fluxy_init.free();
  return;
}
