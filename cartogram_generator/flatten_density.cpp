#include "constants.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include "interpolate_bilinearly.h"
#include <boost/multi_array.hpp>
#include <iostream>
#include <algorithm>
#include <vector>

// Function to calculate the velocity at the grid points (x, y) with x =
// 0.5, 1.5, ..., lx-0.5 and y = 0.5, 1.5, ..., ly-0.5 at time t
void calculate_velocity(double t,
                        FTReal2d &grid_fluxx_init,
                        FTReal2d &grid_fluxy_init,
                        FTReal2d &rho_ft,
                        FTReal2d &rho_init,
                        boost::multi_array<double, 2> *grid_vx,
                        boost::multi_array<double, 2> *grid_vy,
                        const unsigned int lx,
                        const unsigned int ly)
{
  double rho;

#pragma omp parallel for private(rho)
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      rho = rho_ft(0, 0) + (1.0 - t) * (rho_init(i, j) - rho_ft(0,0));
      (*grid_vx)[i][j] = -grid_fluxx_init(i, j) / rho;
      (*grid_vy)[i][j] = -grid_fluxy_init(i, j) / rho;
    }
  }
  return;
}

bool all_points_are_in_domain(double delta_t,
                              boost::multi_array<XYPoint, 2> *proj,
                              boost::multi_array<XYPoint, 2> *v_intp,
                              const unsigned int lx,
                              const unsigned int ly)
{
  // Return false if and only if there exists a point that would be outside
  // [0, lx] x [0, ly]
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      double px = (*proj)[i][j].x;
      double py = (*proj)[i][j].y;
      double vx = (*v_intp)[i][j].x;
      double vy = (*v_intp)[i][j].y;
      if ((px + 0.5*delta_t*vx < 0.0) || (px + 0.5*delta_t*vx > lx)
          || (py + 0.5*delta_t*vy < 0.0) || (py + 0.5*delta_t*vy > ly)) {
        return false;
      }
    }
  }
  return true;
}

// Function to integrate the equations of motion with the fast flow-based
// method
void flatten_density(InsetState *inset_state)
{
  std::cerr << "In flatten_density()" << std::endl;
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  boost::multi_array<XYPoint, 2> &proj = *inset_state->proj();

  // Constants for the numerical integrator
  const double inc_after_acc = 1.1;
  const double dec_after_not_acc = 0.75;
  const double abs_tol = (std::min(lx, ly) * 1e-6);

  // Resize proj multi-array if running for the first time
  if (proj.shape()[0] != lx || proj.shape()[1] != ly) {
    proj.resize(boost::extents[lx][ly]);
  }
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      proj[i][j].x = i + 0.5;
      proj[i][j].y = j + 0.5;
    }
  }
  FTReal2d &rho_ft = *inset_state->ref_to_rho_ft();
  FTReal2d &rho_init = *inset_state->ref_to_rho_init();

  // Allocate memory for the velocity grid
  boost::multi_array<double, 2> grid_vx(boost::extents[lx][ly]);
  boost::multi_array<double, 2> grid_vy(boost::extents[lx][ly]);

  // Prepare Fourier transforms for the flux
  FTReal2d grid_fluxx_init;
  FTReal2d grid_fluxy_init;
  grid_fluxx_init.allocate(lx, ly);
  grid_fluxy_init.allocate(lx, ly);
  grid_fluxx_init.make_fftw_plan(FFTW_RODFT01, FFTW_REDFT01);
  grid_fluxy_init.make_fftw_plan(FFTW_REDFT01, FFTW_RODFT01);

  // eul[i][j] will be the new position of proj[i][j] proposed by a simple
  // Euler step: move a full time interval delta_t with the velocity at time t
  // and position (proj[i][j].x, proj[i][j].y)
  boost::multi_array<XYPoint, 2> eul(boost::extents[lx][ly]);

  // mid[i][j] will be the new displacement proposed by the midpoint
  // method (see comment below for the formula)
  boost::multi_array<XYPoint, 2> mid(boost::extents[lx][ly]);

  // (vx_intp, vy_intp) will be the velocity at position (proj.x, proj.y) at
  // time t
  boost::multi_array<XYPoint, 2> v_intp(boost::extents[lx][ly]);

  // (vx_intp_half, vy_intp_half) will be the velocity at the midpoint
  // (proj.x + 0.5*delta_t*vx_intp, proj.y + 0.5*delta_t*vy_intp) at time
  // t + 0.5*delta_t
  boost::multi_array<XYPoint, 2> v_intp_half(boost::extents[lx][ly]);

  // Initialize the Fourier transforms of gridvx[] and gridvy[] at
  // every point on the lx-times-ly grid at t = 0. We must typecast lx and ly
  // as double-precision numbers. Otherwise the ratios in the denominator
  // will evaluate as zero.
  double dlx = lx;
  double dly = ly;

  // We temporarily insert the Fourier coefficients for the x-components and
  // y-components of the flux vector into grid_fluxx_init and grid_fluxy_init
  for (unsigned int i = 0; i < lx-1; ++i) {
    double di = i;
    for (unsigned int j = 0; j < ly; ++j) {
      double denom = pi * ((di+1)/dlx + (j/(di+1)) * (j/dly) * (dlx/dly));
      grid_fluxx_init(i, j) =
        -rho_ft(i+1, j) / denom;
    }
  }
  for (unsigned int j = 0; j < ly; ++j) {
    grid_fluxx_init(lx-1, j) = 0.0;
  }
  for (unsigned int i=0; i<lx; ++i) {
    double di = i;
    for (unsigned int j = 0; j < ly-1; ++j) {
      double denom = pi * ((di/(j+1)) * (di/dlx) * (dly/dlx) + (j+1)/dly);
      grid_fluxy_init(i, j) =
        -rho_ft(i, j+1) / denom;
    }
  }
  for (unsigned int i=0; i<lx; ++i) {
    grid_fluxy_init(i, ly-1) = 0.0;
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
    calculate_velocity(t,
                       grid_fluxx_init, grid_fluxy_init,
                       rho_ft, rho_init,
                       &grid_vx, &grid_vy,
                       lx, ly);
#pragma omp parallel for
    for (unsigned int i = 0; i < lx; ++i) {
      for (unsigned int j = 0; j < ly; ++j) {

        // We know, either because of the initialization or because of the
        // check at the end of the last iteration, that (proj.x, proj.y)
        // is inside the rectangle [0, lx] x [0, ly]. This fact guarantees
        // that interpolate_bilinearly() is given a point that cannot cause it
        // to fail.
        v_intp[i][j].x =
          interpolate_bilinearly(proj[i][j].x, proj[i][j].y,
                                 &grid_vx, 'x',
                                 lx, ly);
        v_intp[i][j].y =
          interpolate_bilinearly(proj[i][j].x, proj[i][j].y,
                                 &grid_vy, 'y',
                                 lx, ly);
      }
    }
    bool accept = false;
    while (!accept) {

      // Simple Euler step.

#pragma omp parallel for
      for (unsigned int i = 0; i < lx; ++i) {
        for (unsigned int j = 0; j < ly; ++j) {
          eul[i][j].x = proj[i][j].x + v_intp[i][j].x * delta_t;
          eul[i][j].y = proj[i][j].y + v_intp[i][j].y * delta_t;
        }
      }

      // Use "explicit midpoint method"
      // x <- x + delta_t * v_x(x + 0.5*delta_t*v_x(x,y,t),
      //                        y + 0.5*delta_t*v_y(x,y,t),
      //                        t + 0.5*delta_t)
      // and similarly for y.
      calculate_velocity(t + 0.5*delta_t,
                         grid_fluxx_init, grid_fluxy_init,
                         rho_ft, rho_init,
                         &grid_vx, &grid_vy,
                         lx, ly);

      // Make sure we do not pass a point outside [0, lx] x [0, ly] to
      // interpolate_bilinearly(). Otherwise decrease the time step below and
      // try again.
      accept = all_points_are_in_domain(delta_t, &proj, &v_intp, lx, ly);
      if (accept) {

        // Okay, we can run interpolate_bilinearly()

#pragma omp parallel for
        for (unsigned int i = 0; i < lx; ++i) {
          for (unsigned int j = 0; j < ly; ++j) {
            v_intp_half[i][j].x =
              interpolate_bilinearly(
                proj[i][j].x + 0.5*delta_t*v_intp[i][j].x,
                proj[i][j].y + 0.5*delta_t*v_intp[i][j].y,
                &grid_vx, 'x',
                lx, ly);
            v_intp_half[i][j].y =
              interpolate_bilinearly(
                proj[i][j].x + 0.5*delta_t*v_intp[i][j].x,
                proj[i][j].y + 0.5*delta_t*v_intp[i][j].y,
                &grid_vy, 'y',
                lx, ly);
            mid[i][j].x = proj[i][j].x + v_intp_half[i][j].x * delta_t;
            mid[i][j].y = proj[i][j].y + v_intp_half[i][j].y * delta_t;

            // Do not accept the integration step if the maximum squared
            // difference between the Euler and midpoint proposals exceeds
            // abs_tol. Neither should we accept the integration step if one
            // of the positions wandered out of the domain. If one of these
            // problems occurred, decrease the time step.
            if ((mid[i][j].x-eul[i][j].x) * (mid[i][j].x-eul[i][j].x) +
                (mid[i][j].y-eul[i][j].y) * (mid[i][j].y-eul[i][j].y)
                > abs_tol ||
                mid[i][j].x < 0.0 || mid[i][j].x > lx ||
                mid[i][j].y < 0.0 || mid[i][j].y > ly) {
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
    proj = mid;
    delta_t *= inc_after_acc;  // Try a larger step next time
  }
  grid_fluxx_init.destroy_fftw_plan();
  grid_fluxy_init.destroy_fftw_plan();
  grid_fluxx_init.free();
  grid_fluxy_init.free();
  return;
}
