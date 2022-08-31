#ifndef INSET_STATE_H_
#define INSET_STATE_H_

#include "colors.h"
#include "ft_real_2d.h"
#include "geo_div.h"
#include "intersection.h"
#include "xy_point.h"
#include <boost/multi_array.hpp>
#include <cairo/cairo.h>
#include <functional>
#include <map>
#include <nlohmann/json.hpp>
#include <vector>

// TODO: Transfer this struct to colors.h
struct color {
  double r, g, b;
};

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
  std::vector<proj_qd> proj_sequence_;

  Bbox bbox_;  // Bounding box
  fftw_plan bwd_plan_for_rho_{};
  std::unordered_map<std::string, Color> colors_;

  // Unblurred density mean, min, max
  double dens_min_;
  double dens_mean_;
  double dens_max_;

  // Scaling factor to convert albers unit to 512*512 unit
  double latt_const_;

  // Cumulative cartogram projection
  boost::multi_array<XYPoint, 2> cum_proj_;
  fftw_plan fwd_plan_for_rho_{};

  // Geographic divisions in this inset
  std::vector<GeoDiv> geo_divs_;

  // Copy of original data
  std::vector<GeoDiv> geo_divs_original_;

  // Chosen diagonal for each grid cell
  boost::multi_array<int, 2> grid_diagonals_;

  // Map name. Inset position is appended to the name if n_insets > 2.
  std::string inset_name_;
  std::unordered_map<std::string, bool> is_input_target_area_missing_;
  std::unordered_map<std::string, std::string> labels_;
  unsigned int lx_{}, ly_{};  // Lattice dimensions
  unsigned int n_finished_integrations_;
  std::string pos_;  // Position of inset ("C", "T" etc.)
  boost::multi_array<XYPoint, 2> proj_;  // Cartogram projection
  boost::multi_array<XYPoint, 2> identity_proj_;  // Original projection

  // Rasterized density and its Fourier transform
  FTReal2d rho_ft_, rho_init_;
  std::unordered_map<std::string, double> target_areas_;

  // Vertical adjacency graph std::vector<std::vector<intersection> >
  // vertical_adj_;

  // Make default constructor private so that only
  // InsetState(const std::string) can be called as constructor
  InsetState();

public:
  explicit InsetState(std::string);  // Constructor
  void create_delaunay_t();
  void adjust_for_dual_hemisphere();
  void apply_albers_projection();
  void apply_smyth_craster_projection();
  double area_error_at(const std::string &) const;
  void auto_color();  // Automatically color GeoDivs
  Bbox bbox() const;
  void blur_density(
    const double blur_width,
    const bool plot_density,
    const bool is_format_ps);
  void check_topology();
  int chosen_diag(const Point v[4], unsigned int *, bool = false);
  Color color_at(const std::string &) const;
  bool color_found(const std::string &) const;
  bool colors_empty() const;
  unsigned int colors_size() const;
  void create_contiguity_graph(unsigned int);
  void densify_geo_divs();
  void densify_geo_divs_using_delaunay_t();
  void destroy_fftw_plans_for_rho();
  void execute_fftw_bwd_plan() const;
  void execute_fftw_fwd_plan() const;
  void exit_if_not_on_grid_or_edge(Point p1) const;
  void fill_grid_diagonals(bool = false);

  // Density functions
  // Fill map with density, using scanlines
  void fill_with_density(
    const bool plot_density,
    const bool plot_grid_heatmap,
    const bool image_format_ps);

  // Flatten said density with integration
  void flatten_density();
  void flatten_density_with_node_vertices();

  std::vector<GeoDiv> geo_divs() const;
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

  const std::vector<std::vector<intersection>> horizontal_scans(
    unsigned int) const;
  void increment_integration();
  void initialize_cum_proj();
  void initialize_identity_proj();
  void insert_color(const std::string &, Color);
  void insert_color(const std::string &, std::string);
  void insert_label(const std::string &, const std::string &);
  void insert_target_area(const std::string &, double);
  void insert_whether_input_target_area_is_missing(const std::string &, bool);
  std::string inset_name() const;
  nlohmann::json inset_to_geojson(bool, bool = false) const;
  const std::vector<Segment> intersecting_segments(unsigned int) const;
  std::vector<std::vector<intersection>> intersec_with_parallel_to(
    char,
    unsigned int) const;
  bool is_input_target_area_missing(const std::string &) const;
  std::string label_at(const std::string &) const;
  double latt_const() const;
  unsigned int lx() const;
  unsigned int ly() const;
  void make_fftw_plans_for_rho();
  max_area_error_info max_area_error() const;
  std::pair<double, double> max_and_min_grid_cell_area(
    unsigned int cell_width);
  std::pair<Point, Point> max_and_min_grid_cell_area_index(
    unsigned int cell_width);
  unsigned int n_finished_integrations() const;
  unsigned int n_geo_divs() const;
  void normalize_inset_area(double total_cart_target_area, bool = false);
  unsigned long n_points() const;
  unsigned int n_rings() const;
  std::string pos() const;
  void project();
  Point projected_point(Point, bool = false);
  Point projected_point_with_triangulation(Point, bool = false);
  void project_with_cum_proj();
  void project_with_delaunay_t();
  void project_with_triangulation();
  void push_back(const GeoDiv &);
  FTReal2d *ref_to_rho_ft();
  FTReal2d *ref_to_rho_init();
  void remove_tiny_polygons(const double &minimum_polygon_size);
  void replace_target_area(const std::string &, double);
  void rescale_map(unsigned int, bool);
  void revert_smyth_craster_projection();
  void rings_are_simple();
  void set_area_errors();
  void set_grid_dimensions(unsigned int, unsigned int);
  void set_inset_name(const std::string &);
  void set_latt_const(const double);
  void set_pos(const std::string &);
  void simplify(unsigned int);
  void store_original_geo_divs();
  double target_area_at(const std::string &) const;
  bool target_area_is_missing(const std::string &) const;
  double total_inset_area() const;
  double total_target_area() const;
  Polygon transform_to_albers_coor(Polygon edge_points);
  std::array<Point, 3> transformed_triangle(
    const std::array<Point, 3> &,
    bool = false);

  // Apply given function to all points
  void transform_points(const std::function<Point(Point)> &, bool = false);
  void trim_grid_heatmap(cairo_t *cr, double padding);
  std::array<Point, 3> untransformed_triangle(const Point, bool = false);

  // Write all intersections found to an SVG/PS file named
  // "*_intersections_*.svg/ps"
  // void write_intersections_image(unsigned int, const bool);

  // Cairo functions
  void write_cairo_map(const std::string &, bool);
  void write_cairo_polygons_to_png(const std::string &, bool, bool, bool);
  void write_cairo_polygons_to_ps(const std::string &, bool, bool, bool);

  // Functions to write map to eps
  void write_grid_heatmap_data(const std::string filename);

  // Functions to write map to eps
  // void write_density_to_eps(const std::string, const double *);
  void write_grid_heatmap_image(
    const std::string filename,
    const bool plot_equal_area_map,
    const bool image_format_ps,
    const bool crop_polygons);
  void write_grids_to_cairo_surface(cairo_t *cr);
  void write_grid_colors_to_cairo_surface(
    cairo_t *cr,
    bool plot_equal_area_map,
    bool crop_polygons);
  // void write_grid_to_eps(std::ofstream &);
  // void write_intersections_to_eps(unsigned int);
  // void write_map_to_eps(const std::string, const bool);
  void write_map_image(
    const std::string filename,
    const bool fill_polygons,
    const bool plot_grid,
    const bool plot_labels,
    const bool image_format_ps,
    const bool equal_area_map);
  // void write_polygons_to_eps(std::ofstream &, const bool, const bool);
  void write_polygons_to_cairo_surface(
    cairo_t *cr,
    const bool fill_polygons,
    const bool colors,
    const bool plot_equal_area_map);
  void write_labels_to_cairo_surface(cairo_t *cr);
  void write_density_image(
    const std::string filename,
    const double *density,
    const bool plot_grid_heatmap,
    const bool image_format_ps);
  void write_intersections_image(unsigned int res, const bool image_format_ps);
  void write_density_bar_image(
    std::string filename,
    const bool image_format_ps);
  void write_legend_to_cairo_surface(cairo_t *cr, bool equal_area_map);
  void write_quadtree(const std::string &);
  void write_density_to_eps(const std::string &, const double *);
  void write_grid_to_eps(std::ofstream &);
  void write_intersections_to_eps(unsigned int);
  void write_map_to_eps(const std::string &, bool);
  void write_polygons_to_eps(std::ofstream &, bool, bool);
  void write_polygon_points_on_cairo_surface(cairo_t *, color);
};

#endif
