#include "map_state.h"

MapState::MapState(const bool w) : is_world_map_(w)
{
  return;
}

MapState::~MapState()
{
  if (rho_init_.is_allocated()) {
    rho_init_.free_ft();
  }
  if (rho_ft_.is_allocated()) {
    rho_ft_.free_ft();
  }
  fftw_destroy_plan(fwd_plan_);
  fftw_destroy_plan(bwd_plan_);
  return;
}

const unsigned int MapState::n_geo_divs() const
{
  return geo_divs_.size();
}

const std::vector<GeoDiv> MapState::geo_divs() const
{
  return geo_divs_;
}

std::vector<GeoDiv> *MapState::ref_to_geo_divs() {
  return &geo_divs_;
}

const bool MapState::is_world_map() const
{
  return is_world_map_;
}

void MapState::make_grid(const unsigned int x, const unsigned int y)
{
  lx_ = x;
  ly_ = y;
  rho_init_.set_array_size(lx_, ly_);
  rho_init_.allocate_ft();
  rho_ft_.set_array_size(lx_, ly_);
  rho_ft_.allocate_ft();
  fwd_plan_ = fftw_plan_r2r_2d(lx_, ly_,
                               rho_init_.array(), rho_ft_.array(),
                               FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
  bwd_plan_ = fftw_plan_r2r_2d(lx_, ly_,
                               rho_ft_.array(), rho_init_.array(),
                               FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
  return;
}

const unsigned int MapState::lx() const
{
  return lx_;
}

const unsigned int MapState::ly() const
{
  return ly_;
}

FTReal2d *MapState::ref_to_rho_init()
{
  return &rho_init_;
}

FTReal2d *MapState::ref_to_rho_ft()
{
  return &rho_ft_;
}

void MapState::execute_fwd_plan() const
{
  fftw_execute(fwd_plan_);
  return;
}

void MapState::execute_bwd_plan() const
{
  fftw_execute(bwd_plan_);
  return;
}

void MapState::push_back(const GeoDiv gd)
{
  geo_divs_.push_back(gd);
  return;
}
