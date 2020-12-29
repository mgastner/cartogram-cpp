#include "map_state.h"
#include "boost/multi_array.hpp"

// Function to integrate the equations of motion with the fast flow-based
// method.

void flatten_density(MapState *map_state)
{
  std::cout << "In flatten_density()" << std::endl;
  const unsigned int lx = map_state->lx();
  const unsigned int ly = map_state->ly();

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
  boost::multi_array<Point, 2> eul(boost::extents[lx][ly]);

  fftw_destroy_plan(plan_for_grid_fluxx_init);
  fftw_destroy_plan(plan_for_grid_fluxy_init);


  return;
}
