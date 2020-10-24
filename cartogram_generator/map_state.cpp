#include "map_state.h"

MapState::MapState(const bool world_proj) : world(world_proj)
{
  return;
}

MapState::~MapState()
{
  if (rho_init.get_array()) {
    rho_init.ft_free();
  }
  if (rho_ft.get_array()) {
    rho_ft.ft_free();
  }
  return;
}

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

void MapState::make_grid(const unsigned int x, const unsigned int y)
{
  lx = x;
  ly = y;
  rho_init.set_array_size(lx, ly);
  rho_init.ft_alloc();
  rho_ft.set_array_size(lx, ly);
  rho_ft.ft_alloc();
  plan_fwd = fftw_plan_r2r_2d(lx, ly,
                              rho_init.get_array(), rho_ft.get_array(),
                              FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
  plan_bwd = fftw_plan_r2r_2d(lx, ly,
                              rho_ft.get_array(), rho_init.get_array(),
                              FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
  return;
}

unsigned int MapState::get_lx() const
{
  return lx;
}

unsigned int MapState::get_ly() const
{
  return ly;
}

FTReal2d *MapState::ref_to_rho_init()
{
  return &rho_init;
}

FTReal2d *MapState::ref_to_rho_ft()
{
  return &rho_ft;
}

fftw_plan MapState::get_plan_fwd()
{
  return plan_fwd;
}

fftw_plan MapState::get_plan_bwd()
{
  return plan_bwd;
}

void MapState::push_back(const GeoDiv gd)
{
  geo_divs.push_back(gd);
  return;
}
