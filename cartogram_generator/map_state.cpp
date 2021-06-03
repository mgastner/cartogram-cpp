#include "map_state.h"

InsetState::InsetState(std::string inset_pos) :
  inset_pos_(inset_pos),
{
  n_finished_integrations_ = 0;
  fwd_plan_for_rho_ = NULL;
  bwd_plan_for_rho_ = NULL;
  return;
}

InsetState::~InsetState()
{
  if (fwd_plan_for_rho_) {
    fftw_destroy_plan(fwd_plan_for_rho_);
  }
  if (bwd_plan_for_rho_) {
    fftw_destroy_plan(bwd_plan_for_rho_);
  }
  return;
}

unsigned int InsetState::n_geo_divs() const
{
  return geo_divs_.size();
}

const std::vector<GeoDiv> InsetState::geo_divs() const
{
  return geo_divs_;
}

std::vector<GeoDiv> *InsetState::ref_to_geo_divs()
{
  return &geo_divs_;
}

void InsetState::set_geo_divs(std::vector<GeoDiv> geo_divs_new)
{
  geo_divs_.clear();
  geo_divs_ = geo_divs_new;
}

void InsetState::target_areas_insert(const std::string id, const double area)
{
  target_areas.insert(std::pair<std::string, double>(id, area));
  return;
}

void InsetState::colors_insert(const std::string id, std::string color)
{

  // From https://stackoverflow.com/questions/313970/how-to-convert-stdstring-
  // to-lower-case
  std::transform(color.begin(), color.end(), color.begin(), ::tolower);
  Color c(color);
  colors.insert(std::pair<std::string, Color>(id, c));
  return;
}

double InsetState::target_areas_at(const std::string id)
{
  return target_areas.at(id);
}

bool InsetState::target_area_is_missing(const std::string id) const
{

  // We use negative area as indication that GeoDiv has no target area
  return target_areas.at(id) < 0.0;
}

const Color InsetState::colors_at(const std::string id)
{
  return colors.at(id);
}

bool InsetState::colors_empty() const
{
  return colors.empty();
}

void InsetState::make_grid(const unsigned int x, const unsigned int y)
{
  lx_ = x;
  ly_ = y;
  rho_init_.set_array_size(lx_, ly_);
  rho_init_.allocate_ft();
  rho_ft_.set_array_size(lx_, ly_);
  rho_ft_.allocate_ft();
  fwd_plan_for_rho_ =
    fftw_plan_r2r_2d(lx_, ly_,
                     rho_init_.array(), rho_ft_.array(),
                     FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
  bwd_plan_for_rho_ =
    fftw_plan_r2r_2d(lx_, ly_,
                     rho_ft_.array(), rho_init_.array(),
                     FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
  return;
}

unsigned int InsetState::lx() const
{
  return lx_;
}

unsigned int InsetState::ly() const
{
  return ly_;
}

unsigned int InsetState::new_xmin() const
{
  return new_xmin_;
}

unsigned int InsetState::new_ymin() const
{
  return new_ymin_;
}

void InsetState::set_new_xmin(const unsigned int new_xmin)
{
  new_xmin_ = new_xmin;
}

void InsetState::set_new_ymin(const unsigned int new_ymin)
{
  new_ymin_ = new_ymin;
}

double InsetState::map_scale() const
{
  return map_scale_;
}

void InsetState::set_map_scale(const double map_scale)
{
  map_scale_ = map_scale;
}

FTReal2d *InsetState::ref_to_rho_init()
{
  return &rho_init_;
}

FTReal2d *InsetState::ref_to_rho_ft()
{
  return &rho_ft_;
}

void InsetState::execute_fwd_plan() const
{
  fftw_execute(fwd_plan_for_rho_);
  return;
}

void InsetState::execute_bwd_plan() const
{
  fftw_execute(bwd_plan_for_rho_);
  return;
}

void InsetState::push_back(const GeoDiv gd)
{
  geo_divs_.push_back(gd);
  return;
}

unsigned int InsetState::n_finished_integrations() const
{
  return n_finished_integrations_;
}

void InsetState::inc_integration()
{
  n_finished_integrations_ += 1;
}

boost::multi_array<XYPoint, 2> *InsetState::proj()
{
  return &proj_;
}

double InsetState::max_area_err()
{

  // Formula for relative area error:
  // area_on_cartogram / target_area - 1

  double sum_target_area = 0.0;
  double sum_cart_area = 0.0;
  for (auto gd : geo_divs_) {
    if (!target_area_is_missing(gd.id())) {
      sum_target_area += target_areas_at(gd.id());
      sum_cart_area += gd.area();
    }
  }
  double mae = 0.0;
  for (auto gd : geo_divs_) {
    if (!target_area_is_missing(gd.id())) {
      double obj_area =
        target_areas_at(gd.id()) * sum_cart_area / sum_target_area;
      double relative_area_error = gd.area() / obj_area - 1;
      if (relative_area_error < 0) {
        mae = std::max(mae, -relative_area_error);
      } else {
        mae = std::max(mae, relative_area_error);
      }
    }
  }
  std::cout << "max. area err: " << mae << std::endl;
  return mae;
}

void InsetState::set_inset_pos(std::string inset_pos)
{
  inset_pos_ = inset_pos;
}

CartogramInfo::CartogramInfo(std::string v, const bool w, const bool wd2eps) :
  visual_variable_file_(v),
  is_world_map_(w),
  write_density_to_eps_(wd2eps)
{
  return;
}

void CartogramInfo::set_id_header(const std::string id)
{
  id_header_ = id;
  return;
}

const std::string CartogramInfo::id_header() const
{
  return id_header_;
}

const std::string CartogramInfo::visual_variable_file() const
{
  return visual_variable_file_;
}

void CartogramInfo::insert_id_in_visual_variables_file(const std::string id)
{
  ids_in_visual_variables_file_.insert(id);
}

const std::set<std::string> CartogramInfo::ids_in_visual_variables_file() const
{
  return ids_in_visual_variables_file_;
}

bool CartogramInfo::is_world_map() const
{
  return is_world_map_;
}

bool CartogramInfo::trigger_write_density_to_eps() const
{
  return write_density_to_eps_;
}

void CartogramInfo::set_map_name(std::string map_name)
{
  map_name_ = map_name;
}

const std::vector<InsetStates> CartogramInfo::inset_states() const
{
  return inset_states_;
}

void CartogramInfo::push_back(const InsetState is)
{
  inset_states_.push_back(is);
  return;
}
