#include "map_state.h"
#include "geo_div.h"
#include <boost/multi_array.hpp>
#include <iostream>
#include <algorithm>
#include <vector>

#define PI 3.14159265358979323846264338327950288419716939937510
#define INC_AFTER_ACC 1.1
#define DEC_AFTER_NOT_ACC 0.75
#define ABS_TOL (std::min(lx, ly) * 1e-6)

void ffb_calcv (double t,
                FTReal2d grid_fluxx_init,
                FTReal2d grid_fluxy_init,
                FTReal2d rho_ft,
                FTReal2d rho_init,
                boost::multi_array<double, 2> grid_vx,
                boost::multi_array<double, 2> grid_vy,
                const unsigned int lx,
                const unsigned int ly)
{
#pragma omp parallel for
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      double rho = rho_ft(0, 0) + (1.0 - t) * (rho_init(i, j) - rho_ft(0,0));
      grid_vx[i][j] = -grid_fluxx_init(i, j) / rho;
      grid_vy[i][j] = -grid_fluxy_init(i, j) / rho;
    }
  }

  return;
}

// Function to bilinearly interpolate a numerical array grid[0..lx*ly-1]
// whose entries are numbers for the positions:
// x = (0.5, 1.5, ..., lx-0.5), y = (0.5, 1.5, ..., ly-0.5).
// The final argument "zero" can take two possible values: 'x' or 'y'. If
// zero==x, the interpolated function is forced to return 0 if x=0 or x=lx.
// This option is suitable fo interpolating from gridvx because there can be
// no flow through the boundary. If zero==y, the interpolation returns 0 if
// y=0 or y=ly, suitable for gridvy. The unconstrained boundary will be
// determined by continuing the function value at 0.5 (or lx-0.5 or ly-0.5)
// all the way to the edge (i.e. the slope is 0 consistent with a cosine
// transform).
double interpol(double x,
                double y,
                const boost::multi_array<double, 2> grid,
                char zero,
                const unsigned int lx,
                const unsigned int ly)
{
  if (x < 0 || x > lx || y < 0 || y > ly) {
    std::cout << "ERROR: coordinate outside bounding box in interpol()." << "\n";
    std::cout << "x=" << x << ", y=" << y << "\n";
    exit(1);
  }
  if (zero != 'x' && zero != 'y') {
    std::cout << "ERROR: unknown argument zero in interpol()." << "\n";
    exit(1);
  }

  double x0 = std::max(0.0, floor(x + 0.5) - 0.5);
  double x1 = std::min(double(lx), floor(x + 0.5) + 0.5);
  double y0 = std::max(0.0, floor(y + 0.5) - 0.5);
  double y1 = std::min(double(ly), floor(y + 0.5) + 0.5);
  double delta_x = (x - x0) / (x1 - x0);
  double delta_y = (y - y0) / (y1 / y0);

  double fx0y0, fx0y1, fx1y0, fx1y1;
  if ((x<0.5 && y<0.5) || (x<0.5 && zero == 'x') || (y<0.5 && zero == 'y')) {
    fx0y0 = 0.0;
  } else {
    fx0y0 = grid[int(x0)][int(y0)];
  }

  if ((x<0.5 && y>=ly-0.5) || (x<0.5 && zero == 'x') || (y>=ly-0.5 && zero == 'y')) {
    fx0y1 = 0.0;
  } else if (x>=0.5 && y>=ly-0.5 && zero == 'x') {
    fx0y1 = grid[int(x0)][int(ly - 1)];
  } else {
    fx0y1 = grid[int(x0)][int(y1)];
  }

  if ((x>=lx-0.5 && y<0.5) || (x>=lx-0.5 && zero == 'x') || (y<0.5 && zero == 'y')) {
    fx1y0 = 0.0;
  } else if (x>=lx-0.5 && y>=0.5 && zero == 'y') {
    fx1y0 = grid[lx-1][int(y0)];
  } else {
    fx1y0 = grid[int(x1)][int(y0)];
  }

  if ((x>=lx-0.5 && y>=ly-0.5) || (x>=lx-0.5 && zero == 'x') || (y>=ly-0.5 && zero == 'y')) {
    fx1y1 = 0.0;
  } else if (x>=lx-0.5 && y<ly-0.5 && zero == 'y') {
    fx1y1 = grid[lx-1][int(y1)];
  } else if (x<lx-0.5 && y>=ly-0.5 && zero == 'x') {
    fx1y1 = grid[int(x1)][ly - 1];
  } else {
    fx1y1 = grid[int(x1)][int(y1)];
  }

  return (1-delta_x)*(1-delta_y)*fx0y0 + (1-delta_x)*delta_y*fx0y1
         + delta_x*(1-delta_y)*fx1y0 + delta_x*delta_y*fx1y1;
}

// Function to integrate the equations of motion with the fast flow-based
// method.

void flatten_density(MapState *map_state)
{
  std::cout << "In flatten_density()" << std::endl;
  const unsigned int lx = map_state->lx();
  const unsigned int ly = map_state->ly();

  // Define proj and proj2 multi arrays
  boost::multi_array<XYPoint, 2> proj = *map_state->proj();
  boost::multi_array<XYPoint, 2> proj2 = *map_state->proj2();

  // Resize multi arrays if running for the first time
  if (map_state->n_finished_integrations() == 0){
    proj.resize(boost::extents[lx][ly]);
    proj2.resize(boost::extents[lx][ly]);
  }

  for (unsigned int i = 0; i < lx; i++){
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

  // init_griv

  double dlx = lx;
  double dly = ly;

  for (unsigned int i = 0; i < lx; i++) {
    for (unsigned int j = 0; j < ly; j++) {
      rho_ft(i, j) /= 4 * lx * ly;
    }
  }

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

  fftw_execute(plan_for_grid_fluxx_init);
  fftw_execute(plan_for_grid_fluxy_init);

  double t = 0.0;
  double delta_t = 1e-2;
  int iter = 0;

  // Integrate.

  while (t < 1.0) {

    // ffb_calcv(t);
#pragma omp parallel for
    for (unsigned int i = 0; i < lx; i++) {
      for (unsigned int j = 0; j < ly; j++) {
        double rho = rho_ft(0, 0) + (1.0 - t) * (rho_init(i, j) - rho_ft(0,0));
        grid_vx[i][j] = -grid_fluxx_init(i, j) / rho;
        grid_vy[i][j] = -grid_fluxy_init(i, j) / rho;
      }
    }

#pragma omp parallel for
    for (unsigned int i = 0; i < lx; i++) {
      for (unsigned int j = 0; j < ly; j++) {
        v_intp[i][j].x = interpol(proj[i][j].x, proj[i][j].y, grid_vx, 'x', lx, ly);
        v_intp[i][j].y = interpol(proj[i][j].x, proj[i][j].y, grid_vy, 'y', lx, ly);
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

      // ffb_calcv(t + 0.5*delta_t);
      
 #pragma omp parallel for
      for (unsigned int i = 0; i < lx; ++i) {
        for (unsigned int j = 0; j < ly; ++j) {
          double rho = rho_ft(0, 0) + (1.0 - (t)) * (rho_init(i, j) - rho_ft(0,0));
          grid_vx[i][j] = -grid_fluxx_init(i, j) / rho;
          grid_vy[i][j] = -grid_fluxy_init(i, j) / rho;
        }
      }

      // Make sure we do not pass a point outside [0, lx] x [0, ly] to
      // interpol(). Otherwise decrease the time step below and try again.

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

        // OK, we can run interpol().

#pragma omp parallel for
        for (unsigned int i = 0; i < lx; i++) {
          for (unsigned int j = 0; j < ly; j++) {
            v_intp_half[i][j].x = interpol(proj[i][j].x + 0.5*delta_t*v_intp[i][j].x,
                                           proj[i][j].y + 0.5*delta_t*v_intp[i][j].y,
                                           grid_vx, 'x', lx, ly);
            v_intp_half[i][j].y = interpol(proj[i][j].x + 0.5*delta_t*v_intp[i][j].x,
                                           proj[i][j].y + 0.5*delta_t*v_intp[i][j].y,
                                           grid_vy, 'y', lx, ly);
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

  // project() from cartogram.c

  // Calculate displacement vector from proj array
  boost::multi_array<double, 2> xdisp(boost::extents[lx][ly]);
  boost::multi_array<double, 2> ydisp(boost::extents[lx][ly]);
  for (unsigned int i = 0; i < lx; i++){
    for (unsigned int j = 0; j < ly; j++) {
      xdisp[i][j] = proj[i][j].x - i - 0.5;
      ydisp[i][j] = proj[i][j].y - i - 0.5;
    }
  }

  // Insert new positions in a new GeoDiv vector
  std::vector<GeoDiv> new_geo_divs;
  for (auto gd : map_state->geo_divs()){
    GeoDiv new_gd(gd.id());
    for (auto pwh : gd.polygons_with_holes()){

      Polygon ext_ring_old = pwh.outer_boundary();
      CGAL::Polygon_2<Epick> ext_ring;

      for (unsigned int i = 0; i < ext_ring_old.size(); i++){
        ext_ring.push_back(Epick::Point_2(
          interpol(ext_ring_old[i][0], ext_ring_old[i][1], xdisp, 'x', lx, ly) + ext_ring_old[i][0],
          interpol(ext_ring_old[i][0], ext_ring_old[i][1], ydisp, 'y', lx, ly) + ext_ring_old[i][1]
        ));
      }

      std::vector<CGAL::Polygon_2<Epick>> int_ring_v;
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); hci++){
        Polygon old_hole = *hci;
        CGAL::Polygon_2<Epick> int_ring;
        for (unsigned int i = 0; i < old_hole.size(); i++){
          int_ring.push_back(Epick::Point_2(
            interpol(old_hole[i][0], old_hole[i][1], xdisp, 'x', lx, ly) + old_hole[i][0],
            interpol(old_hole[i][0], old_hole[i][1], ydisp, 'y', lx, ly) + old_hole[i][1]
          ));
        }
        int_ring_v.push_back(int_ring);
      }
      const Polygon_with_holes new_pwh(ext_ring, int_ring_v.begin(), int_ring_v.end());
      new_gd.push_back(new_pwh);
    }
    new_geo_divs.push_back(new_gd);
  }

  // Replace old GeoDivs with new ones
  map_state->set_geo_divs(new_geo_divs);


  fftw_destroy_plan(plan_for_grid_fluxx_init);
  fftw_destroy_plan(plan_for_grid_fluxy_init);


  return;
}
