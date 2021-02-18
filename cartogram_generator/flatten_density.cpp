#include "map_state.h"
#include "interpolate_bilinearly.h"
#include <boost/multi_array.hpp>
#include <iostream>
#include <algorithm>
#include <vector>

// Definitions

#define PI 3.14159265358979323846264338327950288419716939937510
#define INC_AFTER_ACC 1.1
#define DEC_AFTER_NOT_ACC 0.75
#define ABS_TOL (std::min(lx, ly) * 1e-6)


// Function to calculate the velocity at the grid points (x, y) with x =
// 0.5, 1.5, ..., lx-0.5 and y = 0.5, 1.5, ..., ly-0.5 at time t.

void ffb_calcv (double t,
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

// Function to integrate the equations of motion with the fast flow-based
// method.

void flatten_density(MapState *map_state)
{
  std::cout << "In flatten_density()" << std::endl;
  const unsigned int lx = map_state->lx();
  const unsigned int ly = map_state->ly();

  // Define proj and proj2 multi arrays
  boost::multi_array<XYPoint, 2> &proj = *map_state->proj();

  // Resize proj multi array if running for the first time
  if (map_state->n_finished_integrations() == 0) {
    proj.resize(boost::extents[lx][ly]);
  }

  for (unsigned int i = 0; i < lx; i++) {
    for (unsigned int j = 0; j < ly; j++) {
      proj[i][j].x = i + 0.5;
      proj[i][j].y = j + 0.5;
    }
  }

  FTReal2d &rho_ft = *map_state->ref_to_rho_ft();
  FTReal2d &rho_init = *map_state->ref_to_rho_init();

  // Allocate memory for the velocity grid
  boost::multi_array<double, 2> grid_vx(boost::extents[lx][ly]);
  boost::multi_array<double, 2> grid_vy(boost::extents[lx][ly]);

  // Prepare Fourier transforms for the flux
  FTReal2d grid_fluxx_init(lx, ly);
  FTReal2d grid_fluxy_init(lx, ly);
  fftw_plan plan_for_grid_fluxx_init =
    fftw_plan_r2r_2d(lx, ly,
                     grid_fluxx_init.array(), grid_fluxx_init.array(),
                     FFTW_RODFT01, FFTW_REDFT01, FFTW_ESTIMATE);
  fftw_plan plan_for_grid_fluxy_init =
    fftw_plan_r2r_2d(lx, ly,
                     grid_fluxy_init.array(), grid_fluxy_init.array(),
                     FFTW_REDFT01, FFTW_RODFT01, FFTW_ESTIMATE);

  // eul[i][j] will be the new position of proj[i][j] proposed by a simple
  // Euler step: move a full time interval delta_t with the velocity at time t
  // and position (proj[i][j].x, proj[i][j].y)
  boost::multi_array<XYPoint, 2> eul(boost::extents[lx][ly]);

  // mid[i*ly+j] will be the new displacement proposed by the midpoint
  // method (see comment below for the formula).
  boost::multi_array<XYPoint, 2> mid(boost::extents[lx][ly]);

  // (vx_intp, vy_intp) will be the velocity at position (proj.x, proj.y) at
  // time t.
  boost::multi_array<XYPoint, 2> v_intp(boost::extents[lx][ly]);

  // (vx_intp_half, vy_intp_half) will be the velocity at the midpoint
  // (proj.x + 0.5*delta_t*vx_intp, proj.y + 0.5*delta_t*vy_intp) at time
  // t + 0.5*delta_t.
  boost::multi_array<XYPoint, 2> v_intp_half(boost::extents[lx][ly]);

  // Initialize the Fourier transforms of gridvx[] and gridvy[] at
  // every point on the lx-times-ly grid at t = 0. After this has
  // finished, we do not need to do any further Fourier transforms for this
  // round of integration

  double dlx = lx;                // We must typecast. Otherwise the ratios in
  double dly = ly;                // the denominator will evaluate as zero.

  // We temporarily insert the Fourier coefficients for the x- and
  // y-components of the flux vector in the arrays grid_fluxx_init[] and
  // grid_fluxy_init[].

  for (unsigned int i = 0; i < lx-1; i++) {
    double di = i;
    for (unsigned int j = 0; j < ly; j++) {
      grid_fluxx_init(i, j) =
        -rho_ft(i+1, j) / (PI * ((di+1)/dlx + (j/(di+1)) * (j/dly) * (dlx/dly)));
    }
  }
  for (unsigned int j = 0; j < ly; j++) {
    grid_fluxx_init(lx-1, j) = 0.0;
  }
  for (unsigned int i=0; i<lx; i++) {
    double di = i;
    for (unsigned int j = 0; j < ly-1; j++) {
      grid_fluxy_init(i, j) =
        -rho_ft(i, j+1) / (PI * ((di/(j+1)) * (di/dlx) * (dly/dlx) + (j+1)/dly));
    }
  }
  for (unsigned int i=0; i<lx; i++) {
    grid_fluxy_init(i, ly-1) = 0.0;
  }

  // Compute the flux vector and store the result in grid_fluxx_init[] and
  // grid_fluxy_init[].

  fftw_execute(plan_for_grid_fluxx_init);
  fftw_execute(plan_for_grid_fluxy_init);

  double t = 0.0;
  double delta_t = 1e-2;                               // Initial time step.
  int iter = 0;

  // Integrate.

  while (t < 1.0) {
    ffb_calcv(t, grid_fluxx_init, grid_fluxy_init,
              rho_ft, rho_init, &grid_vx, &grid_vy, lx, ly);
#pragma omp parallel for
    for (unsigned int i = 0; i < lx; i++) {
      for (unsigned int j = 0; j < ly; j++) {

        // We know, either because of the initialization or because of the
        // check at the end of the last iteration, that (proj.x[k], proj.y[k])
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
      for (unsigned int i = 0; i < lx; i++) {
        for (unsigned int j = 0; j < ly; j++) {
          eul[i][j].x = proj[i][j].x + v_intp[i][j].x * delta_t;
          eul[i][j].y = proj[i][j].y + v_intp[i][j].y * delta_t;
        }
      }

      // Use "explicit midpoint method".
      // x <- x + delta_t * v_x(x + 0.5*delta_t*v_x(x,y,t),
      //                        y + 0.5*delta_t*v_y(x,y,t),
      //                        t + 0.5*delta_t)
      // and similarly for y.

      ffb_calcv(t + 0.5*delta_t, grid_fluxx_init, grid_fluxy_init,
                rho_ft, rho_init, &grid_vx, &grid_vy, lx, ly);

      // Make sure we do not pass a point outside [0, lx] x [0, ly] to
      // interpolate_bilinearly(). Otherwise decrease the time step below and
      // try again.

      accept = true;
      for (unsigned int i = 0; i < lx; i++) {
        for (unsigned int j = 0; j < ly; j++) {
          if (proj[i][j].x + 0.5*delta_t*v_intp[i][j].x < 0.0 ||
                                                          proj[i][j].x + 0.5*delta_t*v_intp[i][j].x > lx ||
              proj[i][j].y + 0.5*delta_t*v_intp[i][j].y < 0.0 ||
                                                          proj[i][j].y + 0.5*delta_t*v_intp[i][j].y > ly) {
            accept = false;
            delta_t *= DEC_AFTER_NOT_ACC;
            break;
          }
        }
        if (!accept) break;
      }

      if (accept) {

        // OK, we can run interpolate_bilinearly().

#pragma omp parallel for
        for (unsigned int i = 0; i < lx; i++) {
          for (unsigned int j = 0; j < ly; j++) {
            v_intp_half[i][j].x =
              interpolate_bilinearly(proj[i][j].x + 0.5*delta_t*v_intp[i][j].x,
                                     proj[i][j].y + 0.5*delta_t*v_intp[i][j].y,
                                     &grid_vx, 'x',
                                     lx, ly);
            v_intp_half[i][j].y =
              interpolate_bilinearly(proj[i][j].x + 0.5*delta_t*v_intp[i][j].x,
                                     proj[i][j].y + 0.5*delta_t*v_intp[i][j].y,
                                     &grid_vy, 'y',
                                     lx, ly);
            mid[i][j].x = proj[i][j].x + v_intp_half[i][j].x * delta_t;
            mid[i][j].y = proj[i][j].y + v_intp_half[i][j].y * delta_t;

            // Do not accept the integration step if the maximum squared
            // difference between the Euler and midpoint proposals exceeds
            // ABS_TOL. Neither should we accept the integration step if one
            // of the positions wandered out of the boundaries. If it
            // happened, decrease the time step.

            if ((mid[i][j].x-eul[i][j].x) * (mid[i][j].x-eul[i][j].x) +
                (mid[i][j].y-eul[i][j].y) * (mid[i][j].y-eul[i][j].y) > ABS_TOL ||
                mid[i][j].x < 0.0 || mid[i][j].x > lx ||
                mid[i][j].y < 0.0 || mid[i][j].y > ly) {
              accept = false;
            }
          }
        }
      }
      if (!accept) {
        delta_t *= DEC_AFTER_NOT_ACC;
      }
    }

    // Control ouput.
    if (iter % 10 == 0) {
      std::cout << "iter = " << iter << ", t = " << t << ", delta_t = " << delta_t << "\n";
    }

    // When we get here, the integration step was accepted.
    t += delta_t;
    iter++;
    boost::multi_array<XYPoint, 2> projtemp = proj;
    proj = mid;
    mid = projtemp;

    delta_t *= INC_AFTER_ACC; // Try a larger step next time.

  }


  fftw_destroy_plan(plan_for_grid_fluxx_init);
  fftw_destroy_plan(plan_for_grid_fluxy_init);


  return;
}
