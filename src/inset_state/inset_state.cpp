#include "inset_state.h"
#include "constants.h"
#include "round_point.h"
#include <cmath>

InsetState::InsetState(std::string pos) : pos_(pos)
{
  n_finished_integrations_ = 0;
  return;
}

void InsetState::create_delaunay_t()
{

  // Get all the points of the map in a unordered_map to remove duplicates
  std::unordered_set<Point> points;

  // Avoid collisions in hash table
  points.reserve(8192);
  points.max_load_factor(0.5);

  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      Polygon ext_ring = pwh.outer_boundary();

      // Get exterior ring coordinates
      for (unsigned int i = 0; i < ext_ring.size(); ++i) {
        points.insert(Point(ext_ring[i][0], ext_ring[i][1]));
      }

      // Get holes of polygon with holes
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
        const Polygon hole = *hci;
        for (unsigned int i = 0; i < hole.size(); ++i) {
          points.insert(Point(hole[i][0], hole[i][1]));
        }
      }
    }
  }

  // Add boundary points of mapping domain
  points.insert(Point(0, 0));
  points.insert(Point(0, ly_));
  points.insert(Point(lx_, 0));
  points.insert(Point(lx_, ly_));

  std::vector<Point> points_vec;

  // Copy points of unordered_set to vector
  std::copy(points.begin(), points.end(), std::back_inserter(points_vec));

  // Create the quadtree and grade it
  Quadtree qt(points_vec, Quadtree::PointMap(), 1);
  int depth = std::max(log2(lx_), log2(ly_));

  std::cerr << "Using Quadtree depth: " << depth << std::endl;

  qt.refine(
    depth,
    9);  // (maximum depth, spltting condition: max number of points per node)
  qt.grade();

  std::cerr << "Quadtree root node bounding box: " << qt.bbox(qt.root())
            << std::endl;

  // Clear corner points from last iteration
  unique_quadtree_corners_.clear();

  // Get unique quadtree corners
  for (Quadtree::Node &node :
       qt.traverse<CGAL::Orthtrees::Leaves_traversal>()) {
        
    // Get bbox of the leaf node
    const Bbox bbox = qt.bbox(node);

    // check if points are between lx_ and ly_
    if (
      bbox.xmin() < 0 || bbox.xmax() > lx_ || bbox.ymin() < 0 ||
      bbox.ymax() > ly_) {
      continue;
    }

    // Insert the four vertices of the bbox into the corners set
    unique_quadtree_corners_.insert(Point(bbox.xmin(), bbox.ymin()));
    unique_quadtree_corners_.insert(Point(bbox.xmax(), bbox.ymax()));
    unique_quadtree_corners_.insert(Point(bbox.xmin(), bbox.ymax()));
    unique_quadtree_corners_.insert(Point(bbox.xmax(), bbox.ymin()));
  }

  // Add boundary points of mapping domain in case they are omitted due to
  // quadtree structure
  unique_quadtree_corners_.insert(Point(0, 0));
  unique_quadtree_corners_.insert(Point(0, ly_));
  unique_quadtree_corners_.insert(Point(lx_, 0));
  unique_quadtree_corners_.insert(Point(lx_, ly_));

  std::cerr << "Number of unique corners: " << unique_quadtree_corners_.size()
            << std::endl;

  // Create the Delaunay triangulation
  Delaunay dt;
  dt.insert(unique_quadtree_corners_.begin(), unique_quadtree_corners_.end());
  proj_qd_.dt = dt;
  std::cerr << "Number of Delaunay triangles: " << dt.number_of_faces()
            << std::endl;
}

double InsetState::area_error_at(const std::string id) const
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
#pragma omp parallel for reduction(min : inset_xmin, inset_ymin) reduction(max : inset_xmax, inset_ymax)
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto bb = pwh.bbox();
      inset_xmin = std::min(bb.xmin(), inset_xmin);
      inset_ymin = std::min(bb.ymin(), inset_ymin);
      inset_xmax = std::max(bb.xmax(), inset_xmax);
      inset_ymax = std::max(bb.ymax(), inset_ymax);
    }
  }
  Bbox inset_bb(inset_xmin, inset_ymin, inset_xmax, inset_ymax);
  return inset_bb;
}

const Color InsetState::color_at(const std::string id) const
{
  return colors_.at(id);
}

bool InsetState::color_found(const std::string id) const
{
  return colors_.count(id);
}

bool InsetState::colors_empty() const
{
  return colors_.empty();
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

void InsetState::increment_integration()
{
  n_finished_integrations_ += 1;
  return;
}

void InsetState::initialize_cum_proj()
{
  cum_proj_.resize(boost::extents[lx_][ly_]);
#pragma omp parallel for
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      cum_proj_[i][j].x = i + 0.5;
      cum_proj_[i][j].y = j + 0.5;
    }
  }
}

void InsetState::insert_color(const std::string id, const Color c)
{
  if (colors_.count(id)) {
    colors_.erase(id);
  }
  colors_.insert(std::pair<std::string, Color>(id, c));
  return;
}

void InsetState::insert_color(const std::string id, std::string color)
{
  if (colors_.count(id)) {
    colors_.erase(id);
  }

  // From https://stackoverflow.com/questions/313970/how-to-convert-stdstring-
  // to-lower-case
  std::transform(color.begin(), color.end(), color.begin(), ::tolower);
  const Color c(color);
  colors_.insert(std::pair<std::string, Color>(id, c));
  return;
}

void InsetState::insert_label(const std::string id, const std::string label)
{
  labels_.insert(std::pair<std::string, std::string>(id, label));
  return;
}

void InsetState::insert_target_area(const std::string id, const double area)
{
  target_areas_.insert(std::pair<std::string, double>(id, area));
  return;
}

void InsetState::insert_whether_input_target_area_is_missing(
  const std::string id,
  const bool is_missing)
{
  is_input_target_area_missing_.insert(
    std::pair<std::string, bool>(id, is_missing));
  return;
}

const std::string InsetState::inset_name() const
{
  return inset_name_;
}

bool InsetState::is_input_target_area_missing(const std::string id) const
{
  return is_input_target_area_missing_.at(id);
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
  fwd_plan_for_rho_ = fftw_plan_r2r_2d(
    lx_,
    ly_,
    rho_init_.as_1d_array(),
    rho_ft_.as_1d_array(),
    FFTW_REDFT10,
    FFTW_REDFT10,
    FFTW_ESTIMATE);
  bwd_plan_for_rho_ = fftw_plan_r2r_2d(
    lx_,
    ly_,
    rho_ft_.as_1d_array(),
    rho_init_.as_1d_array(),
    FFTW_REDFT01,
    FFTW_REDFT01,
    FFTW_ESTIMATE);
  return;
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

unsigned int InsetState::n_finished_integrations() const
{
  return n_finished_integrations_;
}

unsigned int InsetState::n_geo_divs() const
{
  return geo_divs_.size();
}

unsigned long InsetState::n_points() const
{
  unsigned long n_pts = 0;
  for (const auto &gd : geo_divs_) {
    n_pts += gd.n_points();
  }
  return n_pts;
}

unsigned int InsetState::n_rings() const
{
  unsigned int n_rings = 0;
  for (const auto &gd : geo_divs_) {
    n_rings += gd.n_rings();
  }
  return n_rings;
}

const std::string InsetState::pos() const
{
  return pos_;
}

void InsetState::push_back(const GeoDiv gd)
{
  geo_divs_.push_back(gd);
  return;
}

FTReal2d *InsetState::ref_to_rho_ft()
{
  return &rho_ft_;
}

FTReal2d *InsetState::ref_to_rho_init()
{
  return &rho_init_;
}

void InsetState::replace_target_area(const std::string id, const double area)
{
  target_areas_[id] = area;
  return;
}

void InsetState::set_area_errors()
{
  // Formula for relative area error:
  // area_on_cartogram / target_area - 1
  double sum_target_area = 0.0;
  double sum_cart_area = 0.0;
#pragma omp parallel for reduction(+:sum_target_area, sum_cart_area)
  for (const auto &gd : geo_divs_) {
    sum_target_area += target_area_at(gd.id());
    sum_cart_area += gd.area();
  }
  for (const auto &gd : geo_divs_) {
    const double obj_area =
      target_area_at(gd.id()) * sum_cart_area / sum_target_area;
    area_errors_[gd.id()] = std::abs((gd.area() / obj_area) - 1);
  }
  return;
}

void InsetState::set_grid_dimensions(
  const unsigned int lx,
  const unsigned int ly)
{
  lx_ = lx;
  ly_ = ly;
  return;
}

void InsetState::set_inset_name(const std::string inset_name)
{
  inset_name_ = inset_name;
  return;
}

bool InsetState::target_area_is_missing(const std::string id) const
{
  // We use negative area as indication that GeoDiv has no target area
  return target_areas_.at(id) < 0.0;
}

double InsetState::target_area_at(const std::string id) const
{
  return target_areas_.at(id);
}

double InsetState::total_inset_area() const
{
  double total_inset_area = 0.0;
  for (const auto &gd : geo_divs_) {
    total_inset_area += gd.area();
  }
  return total_inset_area;
}

double InsetState::total_target_area() const
{
  double inset_total_target_area = 0;
  for (const auto &geo_div_target_area : target_areas_) {
    inset_total_target_area += geo_div_target_area.second;
  }
  return inset_total_target_area;
}

std::string InsetState::label_at(const std::string id) const
{
  if (labels_.find(id) == labels_.end()) {
    return "";
  }
  return labels_.at(id);
}

void InsetState::transform_points(std::function<Point(Point)> transform_point)
{
  // Iterate over GeoDivs
#pragma omp parallel for
  for (auto &gd : geo_divs_) {

    // Iterate over Polygon_with_holes
    for (auto &pwh : *gd.ref_to_polygons_with_holes()) {

      // Get outer boundary
      auto &outer_boundary = *(&pwh.outer_boundary());

      // Iterate over outer boundary's coordinates
      for (auto &coords_outer : outer_boundary) {

        // Assign outer boundary's coordinates to transformed coordinates
        coords_outer = transform_point(coords_outer);
      }

      // Iterate over holes
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {

        // Iterate over hole's coordinates
        for (auto &coords_hole : *h) {

          // Assign hole's coordinates to transformed coordinates
          coords_hole = transform_point(coords_hole);
        }
      }
    }
  }
  return;
}

void InsetState::round_points()
{

  // Specialise rounded_point with lx_ and ly_
  std::function<Point(Point)> lambda = [lx = lx_, ly = ly_](Point p1) {
    return rounded_point(p1, lx, ly);
  };

  // Apply `lambda` to all points
  transform_points(lambda);
  return;
}
