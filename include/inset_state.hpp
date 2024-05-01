#ifndef INSET_STATE_HPP_
#define INSET_STATE_HPP_

#include "colors.hpp"
#include "ft_real_2d.hpp"
#include "geo_div.hpp"
#include "intersection.hpp"
#include <boost/multi_array.hpp>
#include <cairo/cairo.h>
#include <nlohmann/json.hpp>

struct max_area_error_info {
  double value;
  std::string geo_div;
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

  // Bounding boxes of Quadtree cells
  std::vector<Bbox> quadtree_bboxes_;

  Bbox bbox_;
  fftw_plan bwd_plan_for_rho_{};
  std::unordered_map<std::string, Color> colors_;

  // Unblurred density mean, min, max
  double dens_min_;
  double dens_mean_;
  double dens_max_;

  // Scaling factor to convert equal-area-projection unit to lx*ly unit.
  double latt_const_;

  // Cumulative cartogram projection
  boost::multi_array<Point, 2> cum_proj_;
  fftw_plan fwd_plan_for_rho_{};

  // Geographic divisions in this inset
  std::vector<GeoDiv> geo_divs_;

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
  unsigned int n_finished_integrations_;
  std::string pos_;  // Position of inset ("C", "T" etc.)
  boost::multi_array<Point, 2> proj_;  // Cartogram projection
  boost::multi_array<Point, 2> identity_proj_;  // Original projection

  // Rasterized density, flux and its Fourier transform
  FTReal2d rho_ft_, rho_init_, grid_fluxx_init_, grid_fluxy_init_;
  std::unordered_map<std::string, double> target_areas_;

  // Area errors
  std::vector<double> max_area_errors_;

  // Make default constructor private so that only
  // InsetState(const std::string) can be called as constructor
  InsetState();

public:
  explicit InsetState(std::string);  // Constructor
  void adjust_for_dual_hemisphere();
  void adjust_grid();
  void apply_albers_projection();
  void apply_smyth_craster_projection();

  // Calculate difference between initial area and current area
  double area_drift() const;
  double area_error_at(const std::string &) const;
  void auto_color();  // Automatically color GeoDivs
  Bbox bbox(bool = false) const;
  void blur_density(double, bool);
  double blur_width() const;
  void check_topology();
  int chosen_diag(const Point v[4], unsigned int &, bool = false) const;
  Color color_at(const std::string &) const;
  bool color_found(const std::string &) const;
  bool colors_empty() const;
  unsigned int colors_size() const;
  void create_contiguity_graph(unsigned int);
  void create_delaunay_t();
  void densify_geo_divs();
  void densify_geo_divs_using_delaunay_t();
  void destroy_fftw_plans_for_flux();
  void destroy_fftw_plans_for_rho();
  void execute_fftw_bwd_plan() const;
  void execute_fftw_fwd_plan() const;
  void execute_fftw_plans_for_flux();
  void exit_if_not_on_grid_or_edge(Point p1) const;
  void fill_grid_diagonals(bool = false);

  // Density functions
  void fill_with_density(bool);  // Fill map with density, using scanlines
  void flatten_density();  // Flatten said density with integration
  void flatten_ellipse_density();
  void flatten_density_with_node_vertices();

  const std::vector<GeoDiv> &geo_divs() const;
  std::vector<std::vector<Color>> grid_cell_colors(unsigned int cell_width);
  Polygon grid_cell_edge_points(
    unsigned int x,
    unsigned int y,
    unsigned int cell_width,
    bool plot_equal_area_map);
  double grid_cell_target_area(
    const unsigned int i,
    const unsigned int j,
    const double total_target_area,
    const double total_inset_area);
  double grid_cell_area(
    unsigned int x,
    unsigned int y,
    unsigned int cell_width);
  double grid_cell_area_km(const unsigned int i, const unsigned int j);
  void holes_inside_polygons();
  double grid_cell_target_area_per_km(
    const unsigned int i,
    const unsigned int j,
    const double total_target_area,
    const double total_inset_area);
  Bbox get_bbox_bar(const double bar_width, const double bar_height);

  std::pair<double, unsigned int> get_km_legend_length();
  std::pair<double, unsigned int> get_visual_variable_legend_length();

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
  std::vector<Segment> intersecting_segments(unsigned int) const;
  std::vector<std::vector<intersection>> intersec_with_parallel_to(
    char,
    unsigned int) const;
  bool is_input_target_area_missing(const std::string &) const;
  std::string label_at(const std::string &) const;
  double latt_const() const;
  unsigned int lx() const;
  unsigned int ly() const;
  void make_fftw_plans_for_flux();
  void make_fftw_plans_for_rho();
  void min_ellipses();
  max_area_error_info max_area_error() const;
  std::pair<double, double> max_and_min_grid_cell_area(
    unsigned int cell_width);
  std::pair<Point, Point> max_and_min_grid_cell_area_index(
    unsigned int cell_width);
  unsigned int n_finished_integrations() const;
  unsigned int n_geo_divs() const;
  unsigned long n_points() const;
  unsigned int n_rings() const;
  void normalize_inset_area(double total_cart_target_area, bool = false);
  void normalize_target_area();
  std::string pos() const;
  void project();
  Point projected_point(const Point &, bool = false) const;
  Point projected_point_with_triangulation(const Point &, bool = false) const;
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
  void rescale_map(unsigned int, bool);
  void revert_smyth_craster_projection();
  void rings_are_simple();
  void set_area_errors();
  void set_grid_dimensions(unsigned int, unsigned int);
  void set_geo_divs(std::vector<GeoDiv> new_geo_divs);
  void set_inset_name(const std::string &);
  void store_initial_area();
  void store_initial_target_area();
  void simplify(unsigned int);
  void store_original_geo_divs();
  double target_area_at(const std::string &) const;
  bool target_area_is_missing(const std::string &) const;
  double total_inset_area() const;
  double total_target_area() const;
  Polygon transform_to_equal_area_projection_coor(Polygon edge_points);
  std::array<Point, 3> transformed_triangle(
    const std::array<Point, 3> &,
    bool = false) const;

  // Apply given function to all points
  void transform_points(const std::function<Point(Point)> &, bool = false);
  std::array<Point, 3> untransformed_triangle(const Point &, bool = false)
    const;
  void trim_grid_heatmap(cairo_t *cr, double padding);

  // Cairo functions
  void write_cairo_map(
    const std::string &,
    bool,
    std::unordered_map<Point, Point> = std::unordered_map<Point, Point>());
  void write_cairo_polygons_to_svg(
    const std::string &,
    bool,
    bool,
    bool,
    std::unordered_map<Point, Point> &);

  void write_delaunay_triangles(const std::string &);
  void write_grid_heatmap_data(const std::string filename);

  void write_grid_heatmap_image(
    const std::string filename,
    const bool plot_equal_area_map,
    const bool crop_polygons);
  void write_grids_on_surface(cairo_t *cr);
  void write_grid_colors_on_surface(
    cairo_t *cr,
    bool plot_equal_area_map,
    bool crop_polygons);
  void write_polygons_on_surface(
    cairo_t *cr,
    const bool fill_polygons,
    const bool colors,
    const bool plot_equal_area_map);
  void write_labels_on_surface(cairo_t *cr);
  void write_density_image(
    const std::string filename,
    const double *density,
    const bool plot_pycnophylactic);
  void write_intersections_image(unsigned int res);
  void write_legend_on_surface(cairo_t *cr, bool equal_area_map);
  void write_polygon_points_on_surface(cairo_t *, Color);
  void write_quadtree(const std::string &);
};

#endif  // INSET_STATE_HPP_
