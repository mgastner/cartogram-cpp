#ifndef INSET_STATE_HPP_
#define INSET_STATE_HPP_

#include "colors.hpp"
#include "constants.hpp"
#include "ft_real_2d.hpp"
#include "geo_div.hpp"
#include "intersection.hpp"
#include "nlohmann/json.hpp"
#include "parse_arguments.hpp"
#include "progress_tracker.hpp"
#include "time_tracker.hpp"
#include <boost/multi_array.hpp>
#include <quadtree.hpp>

struct max_area_error_info {
  double value;
  std::string geo_div;

  // Conversion to std::tuple to enable structured bindings
  // this allows:
  // auto [value, geo_div] = max_area_error();
  operator std::tuple<double, std::string>() const
  {
    return std::make_tuple(value, geo_div);
  }
};

struct proj_qd {  // quadtree-delaunay projection
  Delaunay dt;
  std::unordered_map<Point, Point> triangle_transformation;
};

class InsetState
{
private:
  std::unordered_map<std::string, double> area_errors_;

  std::unordered_set<Point> unique_quadtree_corners_;
  proj_qd proj_qd_;

  // Failed constraints
  std::vector<Segment> failed_constraints_dt_projected_;
  std::vector<Segment> failed_constraints_;
  Delaunay og_dt_;

  // New points
  std::unordered_set<Point> points_from_densification_;
  std::unordered_set<Point> points_before_densification_;

  // Bounding boxes of Quadtree cells
  std::vector<Bbox> quadtree_bboxes_;

  Bbox bbox_;
  fftw_plan bwd_plan_for_rho_{};
  std::unordered_map<std::string, Color> colors_;

  Arguments args_;
  TimeTracker timer;

  // File prefix for writing files for that integration
  std::string file_prefix_;

  // Unblurred density mean, min, max
  double dens_min_, dens_mean_, dens_max_, exterior_density_;

  // Scaling factor to convert equal-area-projection unit to lx*ly unit.
  double latt_const_;

  // Cumulative cartogram projection
  boost::multi_array<Point, 2> cum_proj_;
  fftw_plan fwd_plan_for_rho_{};

  // Geographic divisions in this inset
  std::vector<GeoDiv> geo_divs_;

  // Create a map from GeoDiv ID to index in geo_divs_
  std::map<std::string, size_t> geo_divs_id_to_index_;

  // Copy of original data
  std::vector<GeoDiv> geo_divs_original_;

  // Copy to original data to be transform to keep the same points in the final
  // cartogram
  std::vector<GeoDiv> geo_divs_original_transformed_;

  // Chosen diagonal for each grid cell
  boost::multi_array<int, 2> grid_diagonals_;

  // Variable to store initial inset area before integration
  double initial_area_;

  // Store initial total target area before normalization of target area
  // to later be able to normalize inset area by comparing among the insets
  double initial_target_area_;

  // Map name. Inset position is appended to the name if n_insets > 2.
  std::string inset_name_;
  std::unordered_map<std::string, bool> is_input_target_area_missing_;
  std::unordered_map<std::string, std::string> labels_;
  unsigned int lx_{}, ly_{};  // Lattice dimensions
  unsigned int n_fails_during_flatten_density_;
  unsigned int n_finished_integrations_;
  std::string pos_;  // Position of inset ("C", "T" etc.)
  boost::multi_array<Point, 2> proj_;  // Cartogram projection
  boost::multi_array<Point, 2> identity_proj_;  // Original projection

  // Rasterized density, flux and its Fourier transform
  FTReal2d rho_ft_, rho_init_, grid_fluxx_init_, grid_fluxy_init_;
  std::unordered_map<std::string, double> target_areas_;

  // Area errors
  std::vector<double> max_area_errors_;

  // Whether convergence has been reached
  mutable bool converge_{true};

  // Make default constructor private so that only
  // InsetState(const std::string, Arguments) can be called as constructor
  InsetState();

public:
  explicit InsetState(std::string, Arguments);
  void adjust_for_dual_hemisphere();
  void adjust_grid();
  void apply_albers_projection();
  void apply_smyth_craster_projection();

  // Calculate difference between initial area and current area
  double area_expansion_factor() const;
  double area_error_at(const std::string &) const;
  void auto_color();  // Automatically color GeoDivs
  Bbox bbox(bool = false) const;
  void blur_density();
  double blur_width() const;
  void check_topology() const;
  int chosen_diag(const Point v[4], unsigned int &, bool = false) const;

  void cleanup_after_integration();

  Color color_at(const std::string &) const;
  bool color_found(const std::string &) const;
  size_t colors_size() const;
  bool continue_integrating() const;
  void create_and_refine_quadtree();
  void create_contiguity_graph();
  void create_delaunay_t();
  bool converged() const;
  void densify_geo_divs();
  void densify_geo_divs_using_delaunay_t();
  void destroy_fftw_plans_for_flux();
  void destroy_fftw_plans_for_rho();
  void execute_fftw_bwd_plan() const;
  void execute_fftw_fwd_plan() const;
  void execute_fftw_plans_for_flux();
  void exit_if_not_on_grid_or_edge(Point p1) const;

  // Write CSV of time and max_area_error per integration
  void export_time_report() const;

  void fill_grid_diagonals(bool = false);

  // Density functions
  void fill_with_density();
  void fill_with_density_rays();  // Fill map with density, using scanlines
  void fill_with_density_clip();  // Fill map with density, using clipping
  bool flatten_density();  // Flatten said density with integration
  void flatten_ellipse_density();
  void flatten_density_on_square_grid();
  bool flatten_density_on_node_vertices();  // Bool to check if failed

  const std::vector<GeoDiv> &geo_divs() const;
  const GeoDiv &geo_div_at_id(std::string id) const;
  GeoDiv &geo_div_at_id(std::string id);
  Polygon grid_cell_edge_points(
    unsigned int x,
    unsigned int y,
    unsigned int cell_width = plotted_cell_length,
    bool plot_equal_area_map = false) const;
  void holes_inside_polygons() const;

  void increment_n_fails_during_flatten_density();
  void increment_integration();
  double initial_area() const;
  double initial_target_area() const;
  void initialize_cum_proj();
  void initialize_identity_proj();
  void insert_color(const std::string &, const Color &);
  void insert_color(const std::string &, std::string &);
  void insert_label(const std::string &, const std::string &);
  void insert_target_area(const std::string &, double);
  void insert_whether_input_target_area_is_missing(const std::string &, bool);
  std::string inset_name() const;
  nlohmann::json inset_to_geojson(bool, bool = false) const;

  // Function to go from equal area to cartogram
  void integrate(ProgressTracker &);

  std::vector<Segment> intersecting_segments(unsigned int) const;
  std::vector<std::vector<intersection>> intersec_with_parallel_to(
    char,
    unsigned int) const;
  bool is_input_target_area_missing(const std::string &) const;
  void is_simple(const char *caller_func) const;
  std::string label_at(const std::string &) const;
  double latt_const() const;
  unsigned int lx() const;
  unsigned int ly() const;
  void make_fftw_plans_for_flux();
  void make_fftw_plans_for_rho();
  void min_ellipses();
  max_area_error_info max_area_error() const;
  unsigned int n_finished_integrations() const;
  unsigned int n_fails_during_flatten_density() const;
  unsigned int n_geo_divs() const;
  unsigned long n_points() const;
  unsigned int n_rings() const;
  void normalize_inset_area(
    double total_cart_target_area,
    bool equal_area = false,
    bool normalize_original = false);
  void normalize_target_area();
  std::string pos() const;

  void preprocess();
  void prepare_for_integration();

  void print_time_report() const;

  void project();
  void project_with_bilinear_interpolation();
  Point projected_point(const Point &, bool = false) const;
  Point projected_point_with_triangulation(const Point &, bool = false) const;
  void project_point_set(std::unordered_set<Point> &unprojected);
  void project_with_cum_proj();
  void project_with_delaunay_t(bool);
  void project_with_triangulation();
  void push_back(const GeoDiv &);
  FTReal2d &ref_to_fluxx_init();
  FTReal2d &ref_to_fluxy_init();
  FTReal2d &ref_to_rho_ft();
  FTReal2d &ref_to_rho_init();
  void remove_tiny_polygons(const double &minimum_polygon_size);
  void reset_n_finished_integrations();
  void replace_target_area(const std::string &, double);
  void rescale_map();
  void revert_smyth_craster_projection();
  void set_area_errors();
  void set_grid_dimensions(unsigned int, unsigned int);
  void set_geo_divs(std::vector<GeoDiv> new_geo_divs);
  void set_inset_name(const std::string &);
  void simplify(unsigned int);
  void store_initial_area();
  void store_initial_target_area(const double override = 0.0);
  void store_original_geo_divs();
  template <class MetricFn>
  void store_quadtree_cell_corners(const Quadtree<MetricFn> &qt);
  double target_area_at(const std::string &) const;
  bool target_area_is_missing(const std::string &) const;
  double total_inset_area(bool = false) const;
  double total_target_area() const;
  std::array<Point, 3> transformed_triangle(
    const std::array<Point, 3> &,
    bool = false) const;

  // Apply given function to all points
  void transform_points(const std::function<Point(Point)> &, bool = false);
  void transform_polygons(
    const std::function<Polygon(Polygon)> &,
    bool = false);
  void scale_points(double scale_factor, bool project_original = false);
  void move_points(double dx, double dy, bool project_original = false);
  std::array<Point, 3> untransformed_triangle(const Point &, bool = false)
    const;
  void update_delaunay_t();
  void update_file_prefix();
  void update_gd_ids(const std::map<std::string, std::string> &);

  void write_map(
    const std::string &,
    bool,
    bool equal_area_map = false,
    const std::unordered_map<Point, Vector> =
      std::unordered_map<Point, Vector>()) const;

  void write_delaunay_triangles(const std::string &, const bool);
  void write_grid_heatmap_data(const std::string filename);
  void write_density_image(const std::string filename);
  void write_intersections_image();
  void write_quadtree(const std::string &);
};

#endif  // INSET_STATE_HPP_
