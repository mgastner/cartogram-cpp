#include "map_state.h"
#include <fftw3.h>

MapState::MapState(const bool world_proj) : world(world_proj)
{
  return;
}

int MapState::n_geo_divs(void) const
{
  return geo_divs.size();
}

std::vector<GeoDiv> MapState::get_geo_divs(void) const
{
  return geo_divs;
}

std::vector<GeoDiv> *MapState::ref_to_geo_divs(void) {
  return &geo_divs;
}

bool MapState::is_world_map(void) const
{
  return world;
}

void MapState::set_lx(const int i)
{
  lx = i;
  return;
}

void MapState::set_ly(const int i)
{
  ly = i;
  return;
}

int MapState::get_lx(void)
{
  return lx;
}

int MapState::get_ly(void)
{
  return ly;
}

void MapState::allocate_rho(void)
{
  rho_init = (double*) fftw_malloc(lx * ly * sizeof(double));
  rho_ft = (double*) fftw_malloc(lx * ly * sizeof(double));
  plan_fwd = fftw_plan_r2r_2d(lx, ly,
                              rho_init, rho_ft,
                              FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
  plan_bwd = fftw_plan_r2r_2d(lx, ly,
                              rho_ft, rho_init,
                              FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
  return;
}

double *MapState::get_rho_init(void)
{
  return rho_init;
}

double *MapState::get_rho_ft(void)
{
  return rho_ft;
}

fftw_plan MapState::get_plan_fwd(void)
{
  return plan_fwd;
}

void MapState::free_rho(void)
{
  fftw_free(rho_ft);
  fftw_free(rho_init);
  fftw_destroy_plan(plan_fwd);
  fftw_destroy_plan(plan_bwd);
  return;
}

void MapState::push_back(const GeoDiv gd)
{
  geo_divs.push_back(gd);
  return;
}
