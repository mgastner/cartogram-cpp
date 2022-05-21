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

struct max_area_error_info {
  double value;
  std::string geo_div;
};

class InsetState {
private:
  std::unordered_map<std::string, double> area_errors_;
  Bbox bbox_;  // Bounding box
  fftw_plan bwd_plan_for_rho_;
  std::unordered_map<std::string, Color> colors_;

  // Cumulative cartogram projection
  boost::multi_array<XYPoint, 2> cum_proj_;
  fftw_plan fwd_plan_for_rho_;
  std::vector<GeoDiv> geo_divs_;  // Geographic divisions in this inset

  // Chosen diagonal for each graticule cell
  boost::multi_array<int, 2> graticule_diagonals_;

  // Map name. Inset position is appended to the name if n_insets > 2.
  std::string inset_name_;
  std::unordered_map<std::string, bool> is_input_target_area_missing_;
  std::unordered_map<std::string, std::string> labels_;
  unsigned int lx_, ly_;  // Lattice dimensions
  unsigned int new_xmin_, new_ymin_;  // Map translation vector
  unsigned int n_finished_integrations_;
  std::string pos_;  // Position of inset ("C", "T" etc.)
  boost::multi_array<XYPoint, 2> proj_;  // Cartogram projection

  // Rasterized density and its Fourier transform
  FTReal2d rho_ft_, rho_init_;
  std::unordered_map<std::string, double> target_areas_;

  // Vertical adjacency graph
  std::vector<std::vector<intersection> > vertical_adj_;

  // Create cairo surface
  void write_polygons_to_cairo_surface(cairo_t *, const bool,
                                       const bool, const bool);
  // Make default contructor private so that only
  // InsetState(const std::string) can be called as constructor
  InsetState();

public:
  explicit InsetState(const std::string);  // Constructor
  void adjust_for_dual_hemisphere();
  void apply_albers_projection();
  void apply_smyth_craster_projection();
  double area_error_at(const std::string) const;
  void auto_color();  // Automatically color GeoDivs
  Bbox bbox() const;
  void blur_density(const double, bool);
  void check_topology();
  int chosen_diag(const Point v[4], unsigned int *);
  const Color color_at(const std::string) const;
  bool color_found(const std::string) const;
  bool colors_empty() const;
  unsigned int colors_size() const;
  void create_contiguity_graph(unsigned int);
  void densify_geo_divs();
  void destroy_fftw_plans_for_rho();
  void execute_fftw_bwd_plan() const;
  void execute_fftw_fwd_plan() const;
  void exit_if_not_on_grid_or_edge(const Point p1) const;
  void fill_graticule_diagonals();

  // Density functions
  void fill_with_density(bool);  // Fill map with density, using scanlines
  void flatten_density();  // Flatten said density with integration

  const std::vector<GeoDiv> geo_divs() const;
  void holes_inside_polygons();
  const std::vector<std::vector<intersection> > horizontal_scans(
    unsigned int
  ) const;
  void increment_integration();
  void initialize_cum_proj();
  void insert_color(const std::string, const Color);
  void insert_color(const std::string, const std::string);
  void insert_label(const std::string, const std::string);
  void insert_target_area(const std::string, const double);
  void insert_whether_input_target_area_is_missing(
    const std::string,
    const bool
  );
  const std::string inset_name() const;
  nlohmann::json inset_to_geojson(bool) const;
  const std::vector<Segment> intersecting_segments(unsigned int) const;
  std::vector<std::vector<intersection> > intersec_with_parallel_to(
    char,
    unsigned int
  ) const;
  bool is_input_target_area_missing(const std::string) const;
  std::string label_at(const std::string) const;
  unsigned int lx() const;
  unsigned int ly() const;
  void make_fftw_plans_for_rho();
  struct max_area_error_info max_area_error() const;
  void min_ellipses();
  unsigned int n_finished_integrations() const;
  unsigned int n_geo_divs() const;
  void normalize_inset_area(double total_cart_target_area, bool = false);
  unsigned long n_points() const;
  unsigned int n_rings() const;
  const std::string pos() const;
  void project();
  Point projected_point(const Point);
  Point projected_point_with_triangulation(const Point);
  void project_with_triangulation();
  void push_back(const GeoDiv);
  FTReal2d *ref_to_rho_ft();
  FTReal2d *ref_to_rho_init();
  void replace_target_area(const std::string, const double);
  void rescale_map(unsigned int, bool);
  void revert_smyth_craster_projection();
  void rings_are_simple();
  void set_area_errors();
  void set_grid_dimensions(const unsigned int, const unsigned int);
  void set_inset_name(const std::string);
  void simplify(const unsigned int);
  double target_area_at(const std::string) const;
  bool target_area_is_missing(const std::string) const;
  double total_inset_area() const;
  double total_target_area() const;
  std::array<Point, 3> transformed_triangle(const std::array<Point, 3>);

  // Apply given function to all points
  void transform_points(std::function<Point(Point)>);
  std::array<Point, 3> untransformed_triangle(const Point);

  // Cairo functions
  void write_cairo_map(const std::string, const bool);
  void write_cairo_polygons_to_png(
    const std::string,
    const bool,
    const bool,
    const bool
  );
  void write_cairo_polygons_to_ps(
    const std::string,
    const bool,
    const bool,
    const bool
  );

  // Functions to write map to eps
  void write_density_to_eps(const std::string, const double *);
  void write_graticule_to_eps(std::ofstream &);
  void write_intersections_to_eps(unsigned int);
  void write_map_to_eps(const std::string, const bool);
  void write_polygons_to_eps(std::ofstream &, const bool, const bool);
};

#endif
