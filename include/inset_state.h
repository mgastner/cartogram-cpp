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

  // Cumulative cartogram projection
  boost::multi_array<XYPoint, 2> cum_proj_;
  fftw_plan fwd_plan_for_rho_{};

  // Geographic divisions in this inset
  std::vector<GeoDiv> geo_divs_;

  // Copy of original data
  std::vector<GeoDiv> geo_divs_original_;

  // Chosen diagonal for each graticule cell
  boost::multi_array<int, 2> graticule_diagonals_;

  // Variable to store initial inset area before integration
  double initial_area_;

  // Map name. Inset position is appended to the name if n_insets > 2.
  std::string inset_name_;
  std::unordered_map<std::string, bool> is_input_target_area_missing_;
  std::unordered_map<std::string, std::string> labels_;
  unsigned int lx_{}, ly_{};  // Lattice dimensions
  unsigned int n_finished_integrations_;
  std::string pos_;  // Position of inset ("C", "T" etc.)
  boost::multi_array<XYPoint, 2> proj_;  // Cartogram projection

  // Rasterized density and its Fourier transform
  FTReal2d rho_ft_, rho_init_;
  std::unordered_map<std::string, double> target_areas_;

  // Vertical adjacency graph
  std::vector<std::vector<intersection> > vertical_adj_;

  // Create cairo surface
  void write_polygons_to_cairo_surface(cairo_t *, bool, bool, bool);

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
  void blur_density(double, bool);
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
  void fill_graticule_diagonals(bool = false);

  // Density functions
  void fill_with_density(bool);  // Fill map with density, using scanlines
  void flatten_density();  // Flatten said density with integration
  void flatten_density_with_node_vertices();

  std::vector<GeoDiv> geo_divs() const;
  void holes_inside_polygons();
  void increment_integration();
  void initialize_cum_proj();
  void insert_color(const std::string &, Color);
  void insert_color(const std::string &, std::string);
  void insert_label(const std::string &, const std::string &);
  void insert_target_area(const std::string &, double);
  void insert_whether_input_target_area_is_missing(const std::string &, bool);
  std::string inset_name() const;
  nlohmann::json inset_to_geojson(bool, bool = false) const;
  std::vector<Segment> intersecting_segments(unsigned int) const;
  std::vector<std::vector<intersection> > intersec_with_parallel_to(
    char,
    unsigned int) const;
  bool is_input_target_area_missing(const std::string &) const;
  std::string label_at(const std::string &) const;
  unsigned int lx() const;
  unsigned int ly() const;
  void make_fftw_plans_for_rho();
  struct max_area_error_info max_area_error() const;
  unsigned int n_finished_integrations() const;
  unsigned int n_geo_divs() const;
  unsigned long n_points() const;
  unsigned int n_rings() const;
  void normalize_inset_area(double total_cart_target_area, bool = false);
  void normalize_target_area();
  std::string pos() const;
  void project();
  Point projected_point(Point, bool = false);
  Point projected_point_with_triangulation(Point, bool = false);
  void project_with_cum_proj();
  void project_with_delaunay_t();
  void project_with_triangulation();
  void push_back(const GeoDiv &);

  // Calculate difference between initial area and current area
  double area_drift() const;
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
  void store_initial_area();
  void simplify(unsigned int);
  void store_original_geo_divs();
  double target_area_at(const std::string &) const;
  bool target_area_is_missing(const std::string &) const;
  double total_inset_area() const;
  double total_target_area() const;
  std::array<Point, 3> transformed_triangle(
    const std::array<Point, 3> &,
    bool = false);

  // Apply given function to all points
  void transform_points(const std::function<Point(Point)> &, bool = false);
  std::array<Point, 3> untransformed_triangle(Point, bool = false);

  // Cairo functions
  void write_cairo_map(const std::string &, bool);
  void write_cairo_polygons_to_png(const std::string &, bool, bool, bool);
  void write_cairo_polygons_to_ps(const std::string &, bool, bool, bool);

  // Functions to write map to eps
  void write_quadtree(const std::string &);
  void write_density_to_eps(const std::string &, const double *);
  void write_graticule_to_eps(std::ofstream &);
  void write_intersections_to_eps(unsigned int);
  void write_map_to_eps(const std::string &, bool);
  void write_polygons_to_eps(std::ofstream &, bool, bool);
  void write_polygon_points_on_cairo_surface(cairo_t *, color);
};

#endif
