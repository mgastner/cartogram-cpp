#include "inset_state.hpp"
#include "constants.hpp"
#include "csv.hpp"
#include "quadtree.hpp"
#include "triangulation.hpp"
#include <algorithm>

InsetState::InsetState(std::string pos, Arguments args)
    : args_(args), pos_(pos)
{
  initial_area_ = 0.0;
  n_finished_integrations_ = 0;
  n_fails_during_flatten_density_ = 0;
  dens_min_ = 0.0;
  dens_mean_ = 0.0;
  dens_max_ = 0.0;
  latt_const_ = 0.0;
  initial_target_area_ = 0.0;
}

double InsetState::area_error_at(const std::string &id) const
{
  try {
    return area_errors_.at(id);
  } catch (const std::out_of_range &e) {
    std::cerr << "ERROR: Key '" << id << "' not found in area_errors_. "
              << "Exception: " << e.what() << std::endl;
    // Re-throw, or return a default value
    throw;
  }
}

Bbox InsetState::bbox(bool original_bbox) const
{
  auto &geo_divs = original_bbox ? geo_divs_original_transformed_ : geo_divs_;
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
      1 + log2(std::max(lx(), ly()) / default_long_grid_length)) +
    n_fails_during_flatten_density_;
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

bool InsetState::continue_integrating() const
{

  // Calculate all the necessary information to decide whether to continue
  auto [max_area_err, worst_gd] = max_area_error();

  // A GeoDiv is still above our area error threshold
  bool area_error_above_threshold =
    max_area_err > args_.max_permitted_area_error;

  // Area expansion factor is above our threshold
  // i.e. cartogram has become too big or too small
  double area_drift = area_expansion_factor() - 1.0;
  bool area_expansion_factor_above_threshold =
    std::abs(area_drift) > max_permitted_area_drift;

  // If both the above metrics are above our threshold
  bool has_converged =
    !area_error_above_threshold && !area_expansion_factor_above_threshold;

  // Make sure to not continue endlesslely: cap at max_integrations
  bool within_integration_limit = n_finished_integrations() < max_integrations;

  // We continue if we are within the integration limit and have not converged
  bool continue_integration =
    (within_integration_limit && !has_converged) ||
    (n_finished_integrations_ < args_.min_integrations);

  // Actually hasn't converged, just reached integration limit
  if (!within_integration_limit && !has_converged) {
    converge_ = false;
    std::cerr << "ERROR: Could not converge!" << std::endl;
    if (area_error_above_threshold)
      std::cerr << "Max area error above threshold!" << std::endl;
    if (area_expansion_factor_above_threshold)
      std::cerr << "Area expansion factor above threshold!" << std::endl;
  }

  // Print control output (at end of previous integration)
  std::cerr << "Max. area err: " << max_area_err << ", GeoDiv: " << worst_gd
            << std::endl;
  std::cerr << "Current Area: " << geo_div_at_id(worst_gd).area()
            << ", Target Area: " << target_area_at(worst_gd) << std::endl;
  std::cerr << "Area drift: " << area_drift * 100.0 << "%" << std::endl;

  if (continue_integration) {
    // Print next integration information.
    std::cerr << "\nIntegration number " << n_finished_integrations()
              << std::endl;
    std::cerr << "Dimensions : " << lx_ << " " << ly_ << std::endl;
    std::cerr << "Number of Points: " << n_points() << std::endl;
  }
  return continue_integration;
}

Color InsetState::color_at(const std::string &id) const
{
  try {
    return colors_.at(id);
  } catch (const std::out_of_range &e) {
    std::cerr << "ERROR: Key '" << id << "' not found in colors_. "
              << "Exception: " << e.what() << std::endl;
    // Re-throw, or return a default value
    throw;
  }
}

bool InsetState::color_found(const std::string &id) const
{
  return colors_.count(id);
}

size_t InsetState::colors_size() const
{
  return colors_.size();
}

bool InsetState::converged() const
{
  return converge_;
}

bool InsetState::create_delaunay_t()
{
  timer.start("Delaunay Triangulation");

  if (!triang_.build(&qt_locator_, &proj_data_))
    return false;

  timer.stop("Delaunay Triangulation");
  return true;
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

void InsetState::export_time_report() const
{
  // Open CSV to write, and create empty object to write to
  std::ofstream out_file_csv;
  std::string csv_file_name = inset_name_ + "_time_report.csv";
  out_file_csv.open(csv_file_name);
  if (!out_file_csv) {
    std::cerr << "ERROR writing CSV: failed to open " << csv_file_name
              << ".csv" << std::endl;
  }
  // Each string vectors = one row, starting with column names
  std::vector<std::vector<std::string>> csv_rows(n_finished_integrations_ + 1);

  // Column names
  csv_rows[0].push_back("Integration Number");
  csv_rows[0].push_back("Time (s)");
  csv_rows[0].push_back("Max Area Error");

  for (size_t i = 0; i < n_finished_integrations_; i++) {

    // Integration number
    csv_rows[i + 1].push_back(std::to_string(i));

    // Time taken (in seconds)
    std::string timer_task_name = inset_name_ + "_" + std::to_string(i);
    std::string time_in_seconds = std::to_string(
      static_cast<double>(timer.duration(timer_task_name).count()) / 1000.0);
    csv_rows[i + 1].push_back(time_in_seconds);

    // Max area error for that integration
    std::ostringstream oss;
    oss << std::setprecision(16) << max_area_errors_[i];
    csv_rows[i + 1].push_back(oss.str());
  }

  // Write to CSV object, and close file afterwards
  auto writer = csv::make_csv_writer(out_file_csv);
  for (const auto &row : csv_rows) {
    writer << row;
  }
  out_file_csv.close();
}

const std::vector<GeoDiv> &InsetState::geo_divs() const
{
  return geo_divs_;
}

// Const and non-const version of geo_div_at_id
const GeoDiv &InsetState::geo_div_at_id(std::string id) const
{
  try {
    return geo_divs_[geo_divs_id_to_index_.at(id)];
  } catch (const std::out_of_range &e) {
    std::cerr << "ERROR: Key '" << id
              << "' not found in geo_divs_id_to_index_. "
              << "Exception: " << e.what() << std::endl;
    // Re-throw, or return a default value
    throw;
  }
}
GeoDiv &InsetState::geo_div_at_id(std::string id)
{
  try {
    return geo_divs_[geo_divs_id_to_index_.at(id)];
  } catch (const std::out_of_range &e) {
    std::cerr << "ERROR: Key '" << id
              << "' not found in geo_divs_id_to_index_. "
              << "Exception: " << e.what() << std::endl;
    // Re-throw, or return a default value
    throw;
  }
}

void InsetState::create_and_refine_quadtree()
{
  timer.start("Quadtree");

  // Determine target leaf count as 2^-9 times grid area.
  const size_t target_leaf_count =
    static_cast<size_t>((lx_ * ly_) / args_.quadtree_leaf_count_factor);
  std::cerr << "Quadtree target leaf count (pre-grading): "
            << target_leaf_count << std::endl;

  auto get_rho_diff = [&](uint32_t i, uint32_t j, uint32_t size) {
    double rho_min = std::numeric_limits<double>::infinity();
    double rho_max = -std::numeric_limits<double>::infinity();

    const uint32_t x_end = std::min(i + size, lx_);
    const uint32_t y_end = std::min(j + size, ly_);

    for (uint32_t x = i; x < x_end; ++x) {
      for (uint32_t y = j; y < y_end; ++y) {
        const double rho = ref_to_rho_init()(x, y);
        rho_min = std::min(rho_min, rho);
        rho_max = std::max(rho_max, rho);
      }
    }
    return (rho_max - rho_min);
  };

  Quadtree qt(std::max(lx_, ly_), target_leaf_count, get_rho_diff);
  qt.build();

  const size_t n_leaves_bef_grading = qt.num_leaves();

  qt.grade();

  // Store the bounding boxes of the leaf nodes (updates
  // unique_quadtree_corners_)
  store_quadtree_cell_corners(qt);

  // Only keep minimum data to be able to locate the leaves from arbitrary
  // point
  qt_locator_.build(lx_, ly_, qt.nodes());

  // Build fast indexing so that given the corner value, we can easily get the
  // projected value
  proj_data_.reserve(lx_ + 1, ly_ + 1);
  proj_data_.build_fast_indexing(unique_quadtree_corners_);

  std::cerr << "Quadtree nodes pre-grading: " << n_leaves_bef_grading
            << std::endl;
  std::cerr << "Quadtree nodes post-grading: " << qt.num_leaves() << std::endl;

  timer.stop("Quadtree");
}

template <class QuadtreeImp>
void InsetState::store_quadtree_cell_corners(const QuadtreeImp &qt)
{
  unique_quadtree_corners_.clear();
  quadtree_bboxes_.clear();

  unique_quadtree_corners_.reserve(4 * qt.num_leaves());
  quadtree_bboxes_.reserve(qt.num_leaves());

  for (const auto &leaf : qt.leaves()) {
    const uint32_t xmin = leaf.x;
    const uint32_t ymin = leaf.y;
    const uint32_t xmax = xmin + leaf.size;
    const uint32_t ymax = ymin + leaf.size;

    if ((xmax > lx_) + (ymax > ly_))
      continue;

    Bbox bbox(xmin, ymin, xmax, ymax);
    quadtree_bboxes_.push_back(bbox);

    unique_quadtree_corners_.emplace_back(xmin, ymin);
    unique_quadtree_corners_.emplace_back(xmax, ymin);
    unique_quadtree_corners_.emplace_back(xmin, ymax);
    unique_quadtree_corners_.emplace_back(xmax, ymax);
  }

  quadtree_bboxes_.push_back(Bbox(0, 0, lx_, ly_));
  unique_quadtree_corners_.emplace_back(0, 0);
  unique_quadtree_corners_.emplace_back(0, ly_);
  unique_quadtree_corners_.emplace_back(lx_, 0);
  unique_quadtree_corners_.emplace_back(lx_, ly_);

  // NOTE: Avoiding sort + unique to remove the duplicates, because we want the
  // corners of a quadtree cell to be close to each other as they are likely
  // to be used together during projection for better cache locality.
  std::unordered_set<uint64_t> seen;
  seen.reserve(unique_quadtree_corners_.size() * 2);

  std::vector<QuadtreeCorner> out;
  out.reserve(unique_quadtree_corners_.size());

  for (auto &c : unique_quadtree_corners_) {
    uint64_t key = (uint64_t(c.x()) << 32) | c.y();
    if (seen.insert(key).second) {
      out.push_back(c);
    }
  }

  unique_quadtree_corners_.swap(out);

  std::cerr << "Number of unique quadtree corners: "
            << unique_quadtree_corners_.size() << '\n'
            << "Number of quadtree leaf nodes: " << qt.num_leaves() << '\n';
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
  try {
    return is_input_target_area_missing_.at(id);
  } catch (const std::out_of_range &e) {
    std::cerr << "ERROR: Key '" << id
              << "' not found in is_input_target_area_missing_. "
              << "Exception: " << e.what() << std::endl;
    // Re-throw, or return a default value
    throw;
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
  // Previously used to calculate average area error
  // TODO: max_area_error_info should return more information
  // including average and absolute area error
  // double sum_errors = 0.0;
  // size_t count = 0;
  std::string worst_gd = it->first;
  double value = it->second;
  for (const auto &[gd_id, area_error] : area_errors_) {
    if (area_error > value) {
      value = area_error;
      worst_gd = gd_id;
    }
    // sum_errors += area_error;
    // ++count;
  }
  return {value, worst_gd};
}

unsigned int InsetState::n_finished_integrations() const
{
  return n_finished_integrations_;
}

unsigned int InsetState::n_fails_during_flatten_density() const
{
  return n_fails_during_flatten_density_;
}

unsigned int InsetState::n_geo_divs() const
{
  return static_cast<unsigned int>(geo_divs_.size());
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
  return total_inset_area() / initial_area_;
}

void InsetState::print_time_report() const
{
  timer.print_summary_report();
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

void InsetState::set_area_errors()
{
  // Formula for relative area error:
  // | area_on_cartogram / target_area - 1 |
  // However, we must also either
  // - multiply target area with aef or
  // - divide gd.area() by aef
  // To account for the area drift already introduced.
  // For instance, if the actual cartogram is 10% larger than it was initially,
  // And our GeoDiv is 5% larger than it initially was, it has actually become
  // relatively smaller compared to the total cartogram area. Thus, we must
  // accordingly inflate its target area to account for the area drift.
  double aef = area_expansion_factor();

#pragma omp parallel for default(none) firstprivate(aef) \
  shared(geo_divs_, area_errors_)
  for (const auto &gd : geo_divs_) {
    const double obj_area = target_area_at(gd.id()) * aef;
    area_errors_[gd.id()] = std::abs((gd.area() / obj_area) - 1);
  }
}

void InsetState::adjust_grid()
{
  auto [curr_max_area_error, worst_gd] = max_area_error();
  max_area_errors_.push_back(curr_max_area_error);
  // TODO: Change to a more sophisticated grid adjustment strategy
  // (based on a tolerance of area error)
  if (
    n_finished_integrations_ > 4 &&
    curr_max_area_error >=
      0.999 * max_area_errors_[n_finished_integrations_ - 1]) {

    // Multiply grid size with factor
    std::cerr << "Adjusting grid size." << std::endl;
    if (
      lx_ * default_grid_factor > args_.max_allowed_autoscale_grid_length or
      ly_ * default_grid_factor > args_.max_allowed_autoscale_grid_length) {
      std::cerr << "Cannot increase grid size further. ";
      std::cerr << "Grid size exceeds maximum allowed grid length."
                << std::endl;
      return;
    }
    lx_ *= default_grid_factor;
    ly_ *= default_grid_factor;

    const Transformation scale(CGAL::SCALING, default_grid_factor);
    transform_points(scale);

    initial_area_ *= default_grid_factor * default_grid_factor;
    transform_points(scale, true);

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
  timer.set_name(inset_name);
}

void InsetState::store_initial_area()
{
  initial_area_ = total_inset_area();
}

void InsetState::store_initial_target_area(const double override)
{
  initial_target_area_ =
    almost_equal(override, 0.0) ? total_target_area() : override;
}

bool InsetState::target_area_is_missing(const std::string &id) const
{
  // We use negative area as indication that GeoDiv has no target area
  try {
    return target_areas_.at(id) < 0.0;
  } catch (const std::out_of_range &e) {
    std::cerr << "ERROR: Key '" << id << "' not found in target_areas_. "
              << "Exception: " << e.what() << std::endl;
    // Re-throw, or return a default value
    throw;
  }
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

double InsetState::total_inset_area(bool original_area) const
{
  auto &geo_divs = original_area ? geo_divs_original_transformed_ : geo_divs_;
  double total_inset_area = 0.0;
  for (const auto &gd : geo_divs) {
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

void InsetState::scale_points(double scale_factor, bool project_original)
{
  const Transformation scale(CGAL::SCALING, scale_factor);
  transform_points(scale, project_original);
}

void InsetState::set_geo_divs(std::vector<GeoDiv> new_geo_divs)
{
  geo_divs_ = std::move(new_geo_divs);
}

void InsetState::update_file_prefix()
{
  file_prefix_ = inset_name_ + "_" + std::to_string(n_finished_integrations_);
}

void InsetState::update_gd_ids(
  const std::map<std::string, std::string> &gd_id_map)
{
  for (auto &gd : geo_divs_) {
    try {
      gd.update_id(gd_id_map.at(gd.id()));
    } catch (const std::out_of_range &e) {
      std::cerr << "ERROR: Key '" << gd.id() << "' not found in gd_id_map. "
                << "Exception: " << e.what() << std::endl;
      // Re-throw, or return a default value
      throw;
    }
  }
}

void InsetState::move_points(double dx, double dy, bool project_original)
{
  const Transformation translate(
    CGAL::TRANSLATION,
    CGAL::Vector_2<Scd>(dx, dy));
  transform_points(translate, project_original);
}
