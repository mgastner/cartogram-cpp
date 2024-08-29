#include "inset_state.hpp"
#include "constants.hpp"

InsetState::InsetState(std::string pos) : pos_(std::move(pos))
{
  initial_area_ = 0.0;
  n_finished_integrations_ = 0;
  dens_min_ = 0.0;
  dens_mean_ = 0.0;
  dens_max_ = 0.0;
  latt_const_ = 0.0;
  initial_target_area_ = 0.0;
}

double InsetState::area_error_at(const std::string &id) const
{
  return area_errors_.at(id);
}

Bbox InsetState::bbox(bool original_bbox) const
{
  auto &geo_divs = original_bbox ? geo_divs_original_ : geo_divs_;
  // Find joint bounding box for all "polygons with holes" in this inset
  double inset_xmin = dbl_inf;
  double inset_xmax = -dbl_inf;
  double inset_ymin = dbl_inf;
  double inset_ymax = -dbl_inf;
#pragma omp parallel for default(none) shared(geo_divs) \
  reduction(min : inset_xmin, inset_ymin)               \
  reduction(max : inset_xmax, inset_ymax)
  for (const auto &gd : geo_divs) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto bb = pwh.bbox();
      inset_xmin = std::min(bb.xmin(), inset_xmin);
      inset_ymin = std::min(bb.ymin(), inset_ymin);
      inset_xmax = std::max(bb.xmax(), inset_xmax);
      inset_ymax = std::max(bb.ymax(), inset_ymax);
    }
  }
  return {inset_xmin, inset_ymin, inset_xmax, inset_ymax};
}

double InsetState::blur_width() const
{

  // Blur density to speed up the numerics in flatten_density() below.
  // We slowly reduce the blur width so that the areas can reach their
  // target values. In case of failure during flatten_density, we increase
  // the blur width to avoid the flipped Delaunay triangles.
  // TODO: whenever blur_width hits 0, the maximum area error will start
  //       increasing again and eventually lead to an invalid grid
  //       cell error when projecting with triangulation. Investigate
  //       why. As a temporary fix, we set blur_width to be always
  //       positive, regardless of the number of integrations.
  const unsigned int blur_default_pow =
    static_cast<unsigned int>(
      6 + log2(std::max(lx(), ly()) / default_long_grid_length)) +
    n_fails_during_flatten_density_ * 2;
  double blur_width =
    std::pow(2.0, blur_default_pow - (0.5 * n_finished_integrations_));

  // NOTE: Read TODO above
  // if (inset_state.n_finished_integrations() < max_integrations) {
  //   blur_width =
  //     std::pow(2.0, 5 - int(inset_state.n_finished_integrations()));
  // } else {
  //   blur_width = 0.0;
  // }

  std::cerr << "blur_width = " << blur_width << std::endl;
  return blur_width;
}

void InsetState::check_completion() const
{
  auto [value, geo_div] = max_area_error();
  if (value > max_permitted_area_error) {
    std::cerr << "ERROR: Could not converge, max area error beyond limit ("
              << value << ", " << geo_div << ")" << std::endl;
  }
  double area_expansion_factor_ = area_expansion_factor();
  if (std::abs(area_expansion_factor_ - 1.0) > max_permitted_area_expansion) {
    std::cerr << "ERROR: Area drift beyond limit: "
              << (area_expansion_factor_ - 1.0) * 100.0 << "%" << std::endl;
  }
}

Color InsetState::color_at(const std::string &id) const
{
  return colors_.at(id);
}

bool InsetState::color_found(const std::string &id) const
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

void InsetState::create_delaunay_t()
{
  Delaunay dt;
  dt.insert(unique_quadtree_corners_.begin(), unique_quadtree_corners_.end());
  proj_qd_.dt = dt;
  std::cerr << "Number of Delaunay triangles: " << dt.number_of_faces()
            << std::endl;
}

void InsetState::update_delaunay_t()
{
  // Create the Delauany triangulation from the projected quadtree corners
  std::vector<Point> projected_unique_quadtree_corners;
  for (auto &pt : unique_quadtree_corners_) {
    projected_unique_quadtree_corners.push_back(
      proj_qd_.triangle_transformation.at(pt));
  }

  Delaunay dt_projected;
  dt_projected.insert(
    projected_unique_quadtree_corners.begin(),
    projected_unique_quadtree_corners.end());

  // Reverse map is useful to get the original point of the projected vertex
  std::unordered_map<Point, Point> reverse_triangle_transformation;
  reverse_triangle_transformation.reserve(4 * unique_quadtree_corners_.size());
  for (auto &[key, val] : proj_qd_.triangle_transformation) {
    reverse_triangle_transformation[val] = key;
  }

  // Add the chosen diagonal of the projected Delaunay triangles as constraints
  // to the original Delaunay triangulation
  std::vector<std::pair<Point, Point>> constraints;
  for (Delaunay::Finite_faces_iterator fit = dt_projected.finite_faces_begin();
       fit != dt_projected.finite_faces_end();
       ++fit) {
    Face_handle face = fit;
    const Point p1 = face->vertex(0)->point();
    const Point p2 = face->vertex(1)->point();
    const Point p3 = face->vertex(2)->point();

    const Point p1_orig = reverse_triangle_transformation.at(p1);
    const Point p2_orig = reverse_triangle_transformation.at(p2);
    const Point p3_orig = reverse_triangle_transformation.at(p3);

    // Only pick the edge if it is diagonal
    if (p1_orig.x() != p2_orig.x() && p1_orig.y() != p2_orig.y()) {
      constraints.push_back(std::make_pair(p1_orig, p2_orig));
    }
    if (p2_orig.x() != p3_orig.x() && p2_orig.y() != p3_orig.y()) {
      constraints.push_back(std::make_pair(p2_orig, p3_orig));
    }
    if (p3_orig.x() != p1_orig.x() && p3_orig.y() != p1_orig.y()) {
      constraints.push_back(std::make_pair(p3_orig, p1_orig));
    }
  }

  // Inserting range is faster than inserting one by one
  proj_qd_.dt.insert_constraints(constraints.begin(), constraints.end());
}

void InsetState::destroy_fftw_plans_for_flux()
{
  grid_fluxx_init_.destroy_fftw_plan();
  grid_fluxy_init_.destroy_fftw_plan();
}

void InsetState::destroy_fftw_plans_for_rho()
{
  fftw_destroy_plan(fwd_plan_for_rho_);
  fftw_destroy_plan(bwd_plan_for_rho_);
}

void InsetState::execute_fftw_bwd_plan() const
{
  fftw_execute(bwd_plan_for_rho_);
}

void InsetState::execute_fftw_plans_for_flux()
{
  grid_fluxx_init_.execute_fftw_plan();
  grid_fluxy_init_.execute_fftw_plan();
}

void InsetState::execute_fftw_fwd_plan() const
{
  fftw_execute(fwd_plan_for_rho_);
}

const std::vector<GeoDiv> &InsetState::geo_divs() const
{
  return geo_divs_;
}

void InsetState::create_and_store_quadtree_cell_corners()
{
  std::vector<Point> points;

  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const Polygon &ext_ring = pwh.outer_boundary();

      // Get exterior ring coordinates
      points.insert(
        points.end(),
        ext_ring.vertices_begin(),
        ext_ring.vertices_end());

      // Get holes of polygon with holes
      for (const auto &h : pwh.holes()) {
        points.insert(points.end(), h.vertices_begin(), h.vertices_end());
      }
    }
  }

  // Add boundary points of mapping domain
  points.push_back(Point(0, 0));
  points.push_back(Point(0, ly_));
  points.push_back(Point(lx_, 0));
  points.push_back(Point(lx_, ly_));

  // Remove the duplicates from points
  std::unordered_set<Point> unique_points(points.begin(), points.end());

  std::vector<Point> unique_points_vec(
    unique_points.begin(),
    unique_points.end());

  // Create the quadtree and 'grade' it so that neighboring quadtree leaves
  // differ by a depth that can only be 0 or 1.
  Quadtree qt(unique_points_vec, Quadtree::PointMap(), 1);
  const unsigned int depth =
    static_cast<unsigned int>(std::max(log2(lx_), log2(ly_)));
  std::cerr << "Using Quadtree depth: " << depth << std::endl;

  auto can_split = [&depth, &qt, this](const Quadtree::Node &node) -> bool {
    // if the node depth is greater than depth, do not split
    if (node.depth() >= depth) {
      return false;
    }

    auto bbox = qt.bbox(node);
    double rho_min = 1e9;
    double rho_max = -1e9;

    // get the minimum rho_init of the bbox of the node
    for (unsigned int i = bbox.xmin(); i < bbox.xmax(); ++i) {
      for (unsigned int j = bbox.ymin(); j < bbox.ymax(); ++j) {
        if (i >= this->lx() || j >= this->ly()) {
          continue;
        }
        if (i < 0 || j < 0) {
          continue;
        }
        rho_min = std::min(rho_min, this->ref_to_rho_init()(i, j));
        rho_max = std::max(rho_max, this->ref_to_rho_init()(i, j));
      }
    }
    return rho_max - rho_min >
           (0.001 + pow((1.0 / n_finished_integrations_), 2));
  };
  qt.refine(can_split);
  qt.grade();
  std::cerr << "Quadtree root node bounding box: " << qt.bbox(qt.root())
            << std::endl;

  // Clear corner points from last iteration
  unique_quadtree_corners_.clear();

  // Clear the vector of bounding boxes
  quadtree_bboxes_.clear();

  // Get unique quadtree corners
  for (const auto &node : qt.traverse<CGAL::Orthtrees::Leaves_traversal>()) {

    // Get bounding box of the leaf node
    const Bbox bbox = qt.bbox(node);

    // Store the bounding box
    quadtree_bboxes_.push_back(bbox);

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
}

void InsetState::increment_integration()
{
  n_finished_integrations_ += 1;
}

void InsetState::increment_n_fails_during_flatten_density()
{
  n_fails_during_flatten_density_ += 1;
}

void InsetState::initialize_cum_proj()
{
  cum_proj_.resize(boost::extents[lx_][ly_]);

#pragma omp parallel for default(none)
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      cum_proj_[i][j] = Point(i + 0.5, j + 0.5);
    }
  }
}

void InsetState::initialize_identity_proj()
{
  identity_proj_.resize(boost::extents[lx_][ly_]);
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      identity_proj_[i][j] = Point(i + 0.5, j + 0.5);
    }
  }
}

void InsetState::insert_color(const std::string &id, const Color &c)
{
  if (colors_.count(id)) {
    colors_.erase(id);
  }
  colors_.insert({id, c});
}

void InsetState::insert_color(const std::string &id, std::string &color)
{
  if (colors_.count(id)) {
    colors_.erase(id);
  }

  // From
  // https://stackoverflow.com/questions/313970/how-to-convert-stdstring-to-lower-case
  std::transform(color.begin(), color.end(), color.begin(), ::tolower);
  const Color c(color);
  colors_.insert({id, c});
}

void InsetState::insert_label(const std::string &id, const std::string &label)
{
  labels_.insert({id, label});
}

void InsetState::insert_target_area(const std::string &id, const double area)
{
  target_areas_.insert({id, area});
}

void InsetState::insert_whether_input_target_area_is_missing(
  const std::string &id,
  const bool is_missing)
{
  is_input_target_area_missing_.insert({id, is_missing});
}

std::string InsetState::inset_name() const
{
  return inset_name_;
}

double InsetState::initial_area() const
{
  return initial_area_;
}

double InsetState::initial_target_area() const
{
  return initial_target_area_;
}

bool InsetState::is_input_target_area_missing(const std::string &id) const
{
  return is_input_target_area_missing_.at(id);
}

void InsetState::is_simple() const
{
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      if (!pwh.outer_boundary().is_simple()) {
        std::cerr << "ERROR: Outer boundary is not simple for GeoDiv "
                  << gd.id() << std::endl;
        exit(1);
      }
      for (const auto &h : pwh.holes()) {
        if (!h.is_simple()) {
          std::cerr << gd.id() << std::endl;
          std::cerr << "ERROR: Hole is not simple for GeoDiv " << gd.id()
                    << std::endl;
          exit(1);
        }
      }
    }
  }
}

double InsetState::latt_const() const
{
  return latt_const_;
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
    static_cast<int>(lx_),  // fftw_plan_...() uses signed integers.
    static_cast<int>(ly_),
    rho_init_.as_1d_array(),
    rho_ft_.as_1d_array(),
    FFTW_REDFT10,
    FFTW_REDFT10,
    FFTW_ESTIMATE);
  bwd_plan_for_rho_ = fftw_plan_r2r_2d(
    static_cast<int>(lx_),
    static_cast<int>(ly_),
    rho_ft_.as_1d_array(),
    rho_init_.as_1d_array(),
    FFTW_REDFT01,
    FFTW_REDFT01,
    FFTW_ESTIMATE);
}

void InsetState::make_fftw_plans_for_flux()
{
  grid_fluxx_init_.make_fftw_plan(FFTW_RODFT01, FFTW_REDFT01);
  grid_fluxy_init_.make_fftw_plan(FFTW_REDFT01, FFTW_RODFT01);
}

struct max_area_error_info InsetState::max_area_error() const
{
  auto it = area_errors_.begin();
  std::string worst_gd = it->first;
  double value = it->second;
  for (const auto &[gd_id, area_error] : area_errors_) {
    if (area_error > value) {
      value = area_error;
      worst_gd = gd_id;
    }
  }
  std::cerr << "max. area err: " << value << ", GeoDiv: " << worst_gd
            << std::endl;
  std::cerr << "Current Area: "
            << geo_divs_[geo_divs_id_to_index_.at(worst_gd)].area()
            << ", Target Area: " << target_area_at(worst_gd) << std::endl;
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

void InsetState::normalize_target_area()
{
  double ta = total_target_area();

  // Assign normalized target area to GeoDivs
  for (const auto &gd : geo_divs_) {
    double normalized_target_area =
      (target_area_at(gd.id()) / ta) * initial_area_;
    replace_target_area(gd.id(), normalized_target_area);
  }
}

std::string InsetState::pos() const
{
  return pos_;
}

double InsetState::area_expansion_factor() const
{
  double area_expansion_factor_ = total_inset_area() / initial_area_;
  std::cerr << "Area drift: " << (area_expansion_factor_ - 1.0) * 100.0 << "%"
            << std::endl;
  return area_expansion_factor_;
}

void InsetState::push_back(const GeoDiv &gd)
{
  geo_divs_id_to_index_.insert({gd.id(), geo_divs_.size()});
  geo_divs_.push_back(gd);
}

FTReal2d &InsetState::ref_to_fluxx_init()
{
  return grid_fluxx_init_;
}

FTReal2d &InsetState::ref_to_fluxy_init()
{
  return grid_fluxy_init_;
}

FTReal2d &InsetState::ref_to_rho_ft()
{
  return rho_ft_;
}

FTReal2d &InsetState::ref_to_rho_init()
{
  return rho_init_;
}

void InsetState::remove_tiny_polygons(const double &minimum_polygon_size)
{
  const double threshold = total_inset_area() * minimum_polygon_size;
  std::vector<GeoDiv> geo_divs_cleaned;

  // Iterate over GeoDivs
  for (auto &gd : geo_divs_) {
    GeoDiv gd_cleaned(gd.id());

    // Sort polygons with holes according to area
    gd.sort_pwh_descending_by_area();
    const auto &pwhs = gd.polygons_with_holes();

    // Iterate over Polygon_with_holes
    for (unsigned int i = 0; i < pwhs.size(); ++i) {
      if (i == 0 || pwh_area(pwhs[i]) > threshold) {
        gd_cleaned.push_back(pwhs[i]);
      }
    }
    geo_divs_cleaned.push_back(gd_cleaned);
  }
  geo_divs_ = std::move(geo_divs_cleaned);
}

void InsetState::replace_target_area(const std::string &id, const double area)
{
  target_areas_[id] = area;
}

void InsetState::reset_n_finished_integrations()
{
  n_finished_integrations_ = 0;
}

void InsetState::set_area_errors()
{
  // Formula for relative area error:
  // area_on_cartogram / target_area - 1
  double sum_target_area = 0.0;
  double sum_cart_area = 0.0;

#pragma omp parallel for default(none) \
  reduction(+ : sum_target_area, sum_cart_area)
  for (const auto &gd : geo_divs_) {
    sum_target_area += target_area_at(gd.id());
    sum_cart_area += gd.area();
  }

#pragma omp parallel for default(none) shared(sum_cart_area, sum_target_area)
  for (const auto &gd : geo_divs_) {
    const double obj_area =
      target_area_at(gd.id()) * sum_cart_area / sum_target_area;
    area_errors_[gd.id()] = std::abs((gd.area() / obj_area) - 1);
  }
}

void InsetState::adjust_grid()
{
  unsigned int long_grid_length = std::max(lx_, ly_);
  double curr_max_area_error = max_area_error().value;
  unsigned int grid_factor =
    (long_grid_length > default_long_grid_length) ? 2 : default_grid_factor;
  max_area_errors_.push_back(curr_max_area_error);
  if (
    n_finished_integrations_ >= 2 &&
    curr_max_area_error > max_area_errors_[n_finished_integrations_ - 1] &&
    curr_max_area_error > max_area_errors_[n_finished_integrations_ - 2]) {

    // Multiply grid size with factor
    std::cerr << "Adjusting grid size." << std::endl;
    if (
      lx_ * grid_factor > max_allowed_autoscale_grid_length or
      ly_ * grid_factor > max_allowed_autoscale_grid_length) {
      std::cerr << "Cannot increase grid size further. ";
      std::cerr << "Grid size exceeds maximum allowed grid length."
                << std::endl;
      return;
    }
    lx_ *= grid_factor;
    ly_ *= grid_factor;

    // Scale all map coordinates
    const Transformation scale(CGAL::SCALING, grid_factor);
    for (auto &gd : geo_divs_) {
      for (auto &pwh : gd.ref_to_polygons_with_holes()) {
        auto &ext_ring = pwh.outer_boundary();
        ext_ring = transform(scale, ext_ring);
        for (auto &h : pwh.holes()) {
          h = transform(scale, h);
        }
      }
    }

    for (auto &gd : geo_divs_original_) {
      for (auto &pwh : gd.ref_to_polygons_with_holes()) {
        auto &ext_ring = pwh.outer_boundary();
        ext_ring = transform(scale, ext_ring);
        for (auto &h : pwh.holes()) {
          h = transform(scale, h);
        }
      }
    }

    initial_area_ *= grid_factor * grid_factor;

    for (auto &gd : geo_divs_original_transformed_) {
      for (auto &pwh : gd.ref_to_polygons_with_holes()) {
        auto &ext_ring = pwh.outer_boundary();
        ext_ring = transform(scale, ext_ring);
        for (auto &h : pwh.holes()) {
          h = transform(scale, h);
        }
      }
    }

    normalize_target_area();

    destroy_fftw_plans_for_rho();
    destroy_fftw_plans_for_flux();
    ref_to_rho_init().free();
    ref_to_rho_ft().free();
    ref_to_fluxx_init().free();
    ref_to_fluxy_init().free();

    // Reallocate FFTW plans
    ref_to_rho_init().allocate(lx_, ly_);
    ref_to_rho_ft().allocate(lx_, ly_);
    ref_to_fluxx_init().allocate(lx_, ly_);
    ref_to_fluxy_init().allocate(lx_, ly_);
    make_fftw_plans_for_rho();
    make_fftw_plans_for_flux();
    initialize_identity_proj();
    initialize_cum_proj();
    set_area_errors();

    Bbox bb = bbox();
    std::cerr << "New grid dimensions: " << lx_ << " " << ly_
              << " with bounding box\n\t(" << bb.xmin() << ", " << bb.ymin()
              << ", " << bb.xmax() << ", " << bb.ymax() << ")" << std::endl;
  }
}

void InsetState::set_grid_dimensions(
  const unsigned int lx,
  const unsigned int ly)
{
  lx_ = lx;
  ly_ = ly;
}

void InsetState::set_inset_name(const std::string &inset_name)
{
  inset_name_ = inset_name;
}

void InsetState::store_initial_area()
{
  initial_area_ = total_inset_area();
}

void InsetState::store_initial_target_area()
{
  initial_target_area_ = total_target_area();
}

bool InsetState::target_area_is_missing(const std::string &id) const
{
  // We use negative area as indication that GeoDiv has no target area
  return target_areas_.at(id) < 0.0;
}

double InsetState::target_area_at(const std::string &id) const
{
  try {
    return target_areas_.at(id);
  } catch (const std::out_of_range &e) {
    std::cerr << "ERROR: Key '" << id << "' not found in target_areas_. "
              << "Exception: " << e.what() << std::endl;
    // Re-throw, or return a default value
    throw;
  }
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

std::string InsetState::label_at(const std::string &id) const
{
  if (labels_.find(id) == labels_.end()) {
    return "";
  }
  return labels_.at(id);
}

void InsetState::store_original_geo_divs()
{
  geo_divs_original_ = geo_divs_;
  geo_divs_original_transformed_ = geo_divs_;
}

void InsetState::transform_points(
  const std::function<Point(Point)> &transform_point,
  bool project_original)
{

  auto &geo_divs =
    project_original ? geo_divs_original_transformed_ : geo_divs_;

  // Iterate over GeoDivs
#pragma omp parallel for default(none) shared(transform_point, geo_divs)
  for (auto &gd : geo_divs) {

    // Iterate over Polygon_with_holes
    for (auto &pwh : gd.ref_to_polygons_with_holes()) {

      // Get outer boundary
      auto &outer_boundary = pwh.outer_boundary();

      // Iterate over outer boundary's coordinates
      for (auto &coords_outer : outer_boundary) {

        // Assign outer boundary's coordinates to transformed coordinates
        coords_outer = transform_point(coords_outer);
      }

      // Iterate over holes
      for (auto &h : pwh.holes()) {

        // Iterate over hole's coordinates
        for (auto &coords_hole : h) {

          // Assign hole's coordinates to transformed coordinates
          coords_hole = transform_point(coords_hole);
        }
      }
    }
  }
}

void InsetState::set_geo_divs(std::vector<GeoDiv> new_geo_divs)
{
  geo_divs_ = std::move(new_geo_divs);
}
