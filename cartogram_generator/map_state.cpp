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
  return;
}

double *MapState::get_rho_init(void)
{
  return rho_init;
}

void MapState::free_rho(void)
{
  fftw_free(rho_ft);
  fftw_free(rho_init);
  return;
}

void MapState::push_back(const GeoDiv gd)
{
  geo_divs.push_back(gd);
  return;
}
