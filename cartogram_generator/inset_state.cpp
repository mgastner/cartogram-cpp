#include "constants.h"
#include "inset_state.h"

InsetState::InsetState(std::string pos) : pos_(pos)
{
  n_finished_integrations_ = 0;
  return;
}

double InsetState::area_errors_at(const std::string id) const
{
  return area_errors_.at(id);
}

Bbox InsetState::bbox() const
{
  // Find joint bounding for all polygons with holes in this inset
  double inset_xmin = dbl_inf;
  double inset_xmax = -dbl_inf;
  double inset_ymin = dbl_inf;
  double inset_ymax = -dbl_inf;
  for (GeoDiv gd : geo_divs_) {
    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      Bbox bb = pgnwh.bbox();
      inset_xmin = std::min(bb.xmin(), inset_xmin);
      inset_ymin = std::min(bb.ymin(), inset_ymin);
      inset_xmax = std::max(bb.xmax(), inset_xmax);
      inset_ymax = std::max(bb.ymax(), inset_ymax);
    }
  }
  Bbox inset_bb(inset_xmin, inset_ymin, inset_xmax, inset_ymax);
  return inset_bb;
}

bool InsetState::color_found(const std::string id) const
{
  return colors_.count(id);
}

const Color InsetState::colors_at(const std::string id) const
{
  return colors_.at(id);
}

bool InsetState::colors_empty() const
{
  return colors_.empty();
}

void InsetState::colors_insert(const std::string id, const Color c)
{
  if (colors_.count(id)) {
    colors_.erase(id);
  }
  colors_.insert(std::pair<std::string, Color>(id, c));
  return;
}

void InsetState::colors_insert(const std::string id, std::string color)
{
  if (colors_.count(id)) {
    colors_.erase(id);
  }

  // From https://stackoverflow.com/questions/313970/how-to-convert-stdstring-
  // to-lower-case
  std::transform(color.begin(), color.end(), color.begin(), ::tolower);
  Color c(color);
  colors_.insert(std::pair<std::string, Color>(id, c));
  return;
}

unsigned int InsetState::colors_size() const
{
  return colors_.size();
}

void InsetState::destroy_fftw_plans_for_rho()
{
  fftw_destroy_plan(fwd_plan_for_rho_);
  fftw_destroy_plan(bwd_plan_for_rho_);
  return;
}

void InsetState::execute_fftw_bwd_plan() const
{
  fftw_execute(bwd_plan_for_rho_);
  return;
}

void InsetState::execute_fftw_fwd_plan() const
{
  fftw_execute(fwd_plan_for_rho_);
  return;
}

const std::vector<GeoDiv> InsetState::geo_divs() const
{
  return geo_divs_;
}

const std::vector<std::vector<intersection> > InsetState::horizontal_adj()
const
{
  return horizontal_adj_;
}

void InsetState::increment_integration()
{
  n_finished_integrations_ += 1;
  return;
}

const std::string InsetState::inset_name() const
{
  return inset_name_;
}

bool InsetState::is_target_area_missing(const std::string id) const
{
  return is_target_area_missing_.at(id);
}

void InsetState::is_target_area_missing_insert(const std::string id,
                                               const bool is_missing)
{
  is_target_area_missing_.insert(std::pair<std::string, bool>(id,
                                                              is_missing));
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

void InsetState::make_fftw_plans_for_rho()
{
  fwd_plan_for_rho_ =
    fftw_plan_r2r_2d(lx_, ly_,
                     rho_init_.as_1d_array(), rho_ft_.as_1d_array(),
                     FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
  bwd_plan_for_rho_ =
    fftw_plan_r2r_2d(lx_, ly_,
                     rho_ft_.as_1d_array(), rho_init_.as_1d_array(),
                     FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
  return;
}

double InsetState::map_scale() const
{
  return map_scale_;
}

struct max_area_error_info InsetState::max_area_error() const
{
  double value = -dbl_inf;
  std::string worst_gd = "";
  for (const auto &[gd_id, area_error] : area_errors_) {
    if (area_error > value) {
      value = area_error;
      worst_gd = gd_id;
    }
  }
  return {value, worst_gd};
}

unsigned int InsetState::new_xmin() const
{
  return new_xmin_;
}

unsigned int InsetState::new_ymin() const
{
  return new_ymin_;
}

unsigned int InsetState::n_finished_integrations() const
{
  return n_finished_integrations_;
}

unsigned int InsetState::n_geo_divs() const
{
  return geo_divs_.size();
}

double InsetState::non_missing_target_area() const
{
  double sum_non_missing_target_area = 0;
  for (auto gd : geo_divs_) {
    if (!target_area_is_missing(gd.id())) {
      sum_non_missing_target_area += gd.area();
    }
  }
  return sum_non_missing_target_area;
}

unsigned long InsetState::n_points() const
{
  unsigned long n_pts = 0;
  for (const auto gd : geo_divs_) {
    for (const auto pwh : gd.polygons_with_holes()) {
      const Polygon &ext_ring = pwh.outer_boundary();
      n_pts += ext_ring.size();
      for (auto hi = pwh.holes_begin(); hi != pwh.holes_end(); ++hi) {
        n_pts += hi->size();
      }
    }
  }
  return n_pts;
}

unsigned int InsetState::n_polygons() const
{
  unsigned int n_pgns = 0;
  for (const auto gd: geo_divs_) {
    for (const auto pwh : gd.polygons_with_holes()) {
      n_pgns += pwh.number_of_holes() + 1;  // Add 1 for external ring 
    }
  }
  return n_pgns;
}

const std::string InsetState::pos() const
{
  return pos_;
}

boost::multi_array<XYPoint, 2> *InsetState::proj()
{
  return &proj_;
}

void InsetState::push_back(const GeoDiv gd)
{
  geo_divs_.push_back(gd);
  return;
}

std::vector<GeoDiv> *InsetState::ref_to_geo_divs()
{
  return &geo_divs_;
}

boost::multi_array<int, 2> *InsetState::ref_to_graticule_diagonals()
{
  return &graticule_diagonals_;
}

FTReal2d *InsetState::ref_to_rho_ft()
{
  return &rho_ft_;
}

FTReal2d *InsetState::ref_to_rho_init()
{
  return &rho_init_;
}

void InsetState::set_area_errors()
{
  // Formula for relative area error:
  // area_on_cartogram / target_area - 1
  double sum_target_area = 0.0;
  double sum_cart_area = 0.0;
  for (auto gd : geo_divs_) {
    sum_target_area += target_areas_at(gd.id());
    sum_cart_area += gd.area();
  }
  for (auto gd : geo_divs_) {
    double obj_area =
      target_areas_at(gd.id()) * sum_cart_area / sum_target_area;
    area_errors_[gd.id()] = std::abs((gd.area() / obj_area) - 1);
  }
  return;
}

void InsetState::set_geo_divs(std::vector<GeoDiv> geo_divs_new)
{
  geo_divs_.clear();
  geo_divs_ = geo_divs_new;
  return;
}

void InsetState::set_grid_dimensions(unsigned int lx, unsigned int ly)
{
  lx_ = lx;
  ly_ = ly;
  return;
}

void InsetState::set_horizontal_adj(std::vector<std::vector<intersection> > ha)
{
  horizontal_adj_.clear();
  horizontal_adj_ = ha;
  return;
}

void InsetState::set_inset_name(std::string inset_name)
{
  inset_name_ = inset_name;
  return;
}

void InsetState::set_map_scale(const double map_scale)
{
  map_scale_ = map_scale;
  return;
}

void InsetState::set_pos(std::string pos)
{
  pos_ = pos;
  return;
}

void InsetState::set_vertical_adj(std::vector<std::vector<intersection> > va)
{
  vertical_adj_.clear();
  vertical_adj_ = va;
  return;
}
void InsetState::set_xmin(const unsigned int new_xmin)
{
  new_xmin_ = new_xmin;
  return;
}

void InsetState::set_ymin(const unsigned int new_ymin)
{
  new_ymin_ = new_ymin;
  return;
}

bool InsetState::target_area_is_missing(const std::string id) const
{

  // We use negative area as indication that GeoDiv has no target area
  return target_areas_.at(id) < 0.0;
}

double InsetState::target_areas_at(const std::string id) const
{
  return target_areas_.at(id);
}

void InsetState::target_areas_insert(const std::string id, const double area)
{
  target_areas_.insert(std::pair<std::string, double>(id, area));
  return;
}

void InsetState::target_areas_replace(const std::string id, const double area)
{
  target_areas_[id] = area;
  return;
}

double InsetState::total_target_area() const
{
  double inset_total_target_area = 0;
  for(auto &geo_div_target_area : target_areas_) {
    inset_total_target_area += geo_div_target_area.second;
  }
  return inset_total_target_area;
}

const std::vector<std::vector<intersection> > InsetState::vertical_adj() const
{
  return vertical_adj_;
}
