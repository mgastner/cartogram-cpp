#include "inset_state.h"

InsetState::InsetState(std::string pos) : pos_(pos)
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

double InsetState::area_errors_at(const std::string id) const
{
  return area_errors_.at(id);
}

CGAL::Bbox_2 InsetState::bbox() const
{
  double inset_xmin = 180;
  double inset_ymin = 90;
  double inset_xmax = -180;
  double inset_ymax = -90;
  for (GeoDiv gd : geo_divs_) {
    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      CGAL::Bbox_2 pgnwh_bbox = pgnwh.bbox();
      inset_xmin = std::min(pgnwh_bbox.xmin(), inset_xmin);
      inset_ymin = std::min(pgnwh_bbox.ymin(), inset_ymin);
      inset_xmax = std::max(pgnwh_bbox.xmax(), inset_xmax);
      inset_ymax = std::max(pgnwh_bbox.ymax(), inset_ymax);
    }
  }
  CGAL::Bbox_2 inset_bbox(inset_xmin, inset_ymin, inset_xmax, inset_ymax);
  return inset_bbox;
}

bool InsetState::color_found(const std::string id) const
{
  return colors_.count(id);
}

const Color InsetState::colors_at(const std::string id)
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

double InsetState::inset_total_target_area()
{
  double inset_total_target_area = 0;
  for(auto &geo_div_target_area : target_areas_) {
    inset_total_target_area += geo_div_target_area.second;
  }
  return inset_total_target_area;
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

void InsetState::execute_bwd_plan() const
{
  fftw_execute(bwd_plan_for_rho_);
  return;
}

void InsetState::execute_fwd_plan() const
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

unsigned int InsetState::lx() const
{
  return lx_;
}

unsigned int InsetState::ly() const
{
  return ly_;
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

double InsetState::map_scale() const
{
  return map_scale_;
}

double InsetState::max_area_error() const
{
  double mae = 0.0;

  for (auto const& [gd_id, area_error] : area_errors_) {
    mae = std::max(mae, area_error);
  }
  std::cout << "max. area err: " << mae << std::endl << std::endl;
  return mae;
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
      area_errors_[gd.id()] = relative_area_error;
    }
  }
  return;
}

double InsetState::cart_area()
{
  double sum_cart_area = 0;
  for (auto gd : geo_divs_) {
    if (!target_area_is_missing(gd.id())) {
      sum_cart_area += gd.area();
    }
  }
  return sum_cart_area;
}

void InsetState::set_geo_divs(std::vector<GeoDiv> geo_divs_new)
{
  geo_divs_.clear();
  geo_divs_ = geo_divs_new;
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

double InsetState::target_areas_at(const std::string id)
{
  return target_areas_.at(id);
}

void InsetState::target_areas_insert(const std::string id, const double area)
{
  target_areas_.insert(std::pair<std::string, double>(id, area));
  return;
}

const std::vector<std::vector<intersection> > InsetState::vertical_adj() const
{
  return vertical_adj_;
}

void InsetState::calculate_bbox() 
{
  GeoDiv gd0 = geo_divs()[0];
  std::vector<Polygon_with_holes> pwhs = gd0.polygons_with_holes();
  CGAL::Bbox_2 bb0 = pwhs[0].bbox();
  double map_xmin = bb0.xmin();
  double map_xmax = bb0.xmax();
  double map_ymin = bb0.ymin();
  double map_ymax = bb0.ymax();

  // Expand bounding box to enclose all GeoDivs
  for (auto gd : geo_divs()) {
    for (auto pwh : gd.polygons_with_holes()) {
      CGAL::Bbox_2 bb = pwh.bbox();
      map_xmin = std::min(map_xmin, bb.xmin());
      map_ymin = std::min(map_ymin, bb.ymin());
      map_xmax = std::max(map_xmax, bb.xmax());
      map_ymax = std::max(map_ymax, bb.ymax());
    }
  }

  bbox_ = {map_xmin, map_ymin, map_xmax, map_ymax};
}

CGAL::Bbox_2 InsetState::bbox()
{
  return bbox_;
}