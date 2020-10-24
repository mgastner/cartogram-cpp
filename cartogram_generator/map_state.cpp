#include "map_state.h"

MapState::MapState(const bool world_proj) : world(world_proj)
{
  return;
}

// MapState::~MapState()
// {
//   if (rho_init) {
//     fftw_free(rho_init);
//   }
//   if (rho_ft) {
//     fftw_free(rho_ft);
//   }
//   return;
// }

unsigned int MapState::n_geo_divs() const
{
  return geo_divs.size();
}

std::vector<GeoDiv> MapState::get_geo_divs() const
{
  return geo_divs;
}

std::vector<GeoDiv> *MapState::ref_to_geo_divs() {
  return &geo_divs;
}

bool MapState::is_world_map() const
{
  return world;
}

// void MapState::make_grid(const unsigned int x, const unsigned int y)
// {
//   lx = x;
//   ly = y;
//   rho_init = (double*) fftw_malloc(lx * ly * sizeof(double));
//   rho_ft = (double*) fftw_malloc(lx * ly * sizeof(double));
//   plan_fwd = fftw_plan_r2r_2d(lx, ly,
//                               rho_init, rho_ft,
//                               FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
//   plan_bwd = fftw_plan_r2r_2d(lx, ly,
//                               rho_ft, rho_init,
//                               FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
//   return;
// }

unsigned int MapState::get_lx() const
{
  return lx;
}

unsigned int MapState::get_ly() const
{
  return ly;
}

// double *MapState::get_rho_init()
// {
//   return rho_init;
// }

// double *MapState::get_rho_ft()
// {
//   return rho_ft;
// }

fftw_plan MapState::get_plan_fwd()
{
  return plan_fwd;
}

void MapState::push_back(const GeoDiv gd)
{
  geo_divs.push_back(gd);
  return;
}
