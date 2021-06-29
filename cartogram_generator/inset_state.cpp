#include "inset_state.h"

InsetState::InsetState(std::string pos) :
  pos_(pos)
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

  if (colors.count(id)) {
    colors.erase(id);
  }

  // From https://stackoverflow.com/questions/313970/how-to-convert-stdstring-
  // to-lower-case
  std::transform(color.begin(), color.end(), color.begin(), ::tolower);
  Color c(color);
  colors.insert(std::pair<std::string, Color>(id, c));
  return;
}

void InsetState::colors_insert(const std::string id, const Color c)
{

  if (colors.count(id)) {
    colors.erase(id);
  }
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

bool InsetState::color_found(const std::string id) const
{
  return colors.count(id);
}

unsigned int InsetState::colors_size() const
{
  return colors.size();
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

void InsetState::set_area_errs()
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

  for (auto gd : geo_divs_) {
    if (!target_area_is_missing(gd.id())) {
      double obj_area =
        target_areas_at(gd.id()) * sum_cart_area / sum_target_area;
      double relative_area_error = std::abs( (gd.area() / obj_area) - 1);
      area_errs[gd.id()] = relative_area_error;
    }
  }
}

double InsetState::area_errs_at(const std::string id) const
{
  return area_errs.at(id);
}

double InsetState::max_area_err() const
{
  double mae = 0.0;

  for (auto const& [gd_id, area_err] : area_errs) {
    mae = std::max(mae, area_err);
  }

  std::cout << "max. area err: " << mae << std::endl << std::endl;
  return mae;
}

void InsetState::set_pos(std::string pos)
{
  pos_ = pos;
}

const std::string InsetState::pos() const
{
  return pos_;
}

void InsetState::set_inset_name(std::string inset_name)
{
  inset_name_ = inset_name;
}

const std::string InsetState::inset_name() const
{
  return inset_name_;
}

void InsetState::set_horizontal_adj(std::vector<std::vector<intersection> > ha)
{
  horizontal_adj_.clear();
  horizontal_adj_ = ha;
}

void InsetState::set_vertical_adj(std::vector<std::vector<intersection> > va)
{
  vertical_adj_.clear();
  vertical_adj_ = va;
}

const std::vector<std::vector<intersection> > InsetState::horizontal_adj() const
{
  return horizontal_adj_;
}

const std::vector<std::vector<intersection> > InsetState::vertical_adj() const
{
  return vertical_adj_;
}
