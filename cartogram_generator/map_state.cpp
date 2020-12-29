#include "map_state.h"

MapState::MapState(std::string v, const bool w) :
  visual_variable_file_(v),
  is_world_map_(w)
{
  return;
}

MapState::~MapState()
{
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

std::vector<GeoDiv> *MapState::ref_to_geo_divs()
{
  return &geo_divs_;
}

void MapState::target_areas_insert(const std::string id, const double area)
{
  target_areas.insert(std::pair<std::string, double>(id, area));
  return;
}

void MapState::colors_insert(const std::string id, std::string color)
{

  // from https://stackoverflow.com/questions/313970/how-to-convert-stdstring-to-lower-case
  std::transform(color.begin(), color.end(), color.begin(), ::tolower);
  Color c(color);
  colors.insert(std::pair<std::string, Color>(id, c));
  return;
}

const double MapState::target_areas_at(const std::string id)
{
  return target_areas.at(id);
}

const Color MapState::colors_at(const std::string id)
{
  return colors.at(id);
}
const bool MapState::colors_empty() const
{
  return colors.empty();
}

void MapState::set_id_header(const std::string id)
{
  id_header_ = id;
  return;
}

const std::string MapState::id_header() const
{
  return id_header_;
}

const std::string MapState::visual_variable_file() const
{
  return visual_variable_file_;
}

void MapState::insert_id_in_visual_variables_file(const std::string id)
{
  ids_in_visual_variables_file_.insert(id);
}

const std::set<std::string> MapState::ids_in_visual_variables_file() const
{
  return ids_in_visual_variables_file_;
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
