#ifndef INSET_STATE_H_
#define INSET_STATE_H_

#include "colors.h"
#include "ft_real_2d.h"
#include "geo_div.h"
#include "xy_point.h"
#include "intersection.h"
#include <vector>
#include <boost/multi_array.hpp>
#include <map>

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
  double map_scale_;  // Lattice spacing in terms of equal-area coordinates
  unsigned int new_xmin_, new_ymin_;  // Map translation vector
  unsigned int n_finished_integrations_;
  std::string pos_;  // Position of inset ("C", "T" etc.)
  boost::multi_array<XYPoint, 2> proj_;  // Cartogram projection

  // Rasterized density and its Fourier transform
  FTReal2d rho_init_, rho_ft_;
  std::unordered_map<std::string, double> target_areas_;

  // Vertical adjacency graph
  std::vector<std::vector<intersection> > vertical_adj_;

  // Make default contructor private so that only
  // InsetState(const std::string) can be called as constructor
  InsetState();

public:
  explicit InsetState(const std::string);  // Constructor
  double area_errors_at(const std::string) const;
  void auto_color();  // Automatically color GeoDivs
  Bbox bbox() const;
  const Color colors_at(const std::string) const;
  bool colors_empty() const;
  bool color_found(const std::string id) const;
  void colors_insert(const std::string, const Color);
  void colors_insert(const std::string, const std::string);
  unsigned int colors_size() const;
  void create_contiguity_graph(unsigned int);
  void densify_geo_divs();
  void destroy_fftw_plans_for_rho();
  void execute_fftw_bwd_plan() const;
  void execute_fftw_fwd_plan() const;
  void fill_with_density(bool);  // Fill map with density, using scanlines
  const std::vector<GeoDiv> geo_divs() const;
  const std::vector<std::vector<intersection> >
    horizontal_scans(unsigned int) const;
  void increment_integration();
  void initialize_cum_proj();
  const std::string inset_name() const;
  const std::vector<Segment> intersecting_segments(unsigned int) const;
  bool is_input_target_area_missing(const std::string) const;
  void is_input_target_area_missing_insert(const std::string, const bool);
  std::string labels_at(const std::string) const;
  void labels_insert(const std::string, const std::string);
  unsigned int lx() const;
  unsigned int ly() const;
  void make_fftw_plans_for_rho();
  double map_scale() const;
  struct max_area_error_info max_area_error() const;
  double max_area_error_stdcerr() const;
  unsigned int new_xmin() const;
  unsigned int new_ymin() const;
  unsigned int n_finished_integrations() const;
  unsigned int n_geo_divs() const;
  unsigned long n_points() const;
  unsigned int n_rings() const;
  const std::string pos() const;
  void push_back(const GeoDiv);
  boost::multi_array<XYPoint, 2> *ref_to_cum_proj();
  std::vector<GeoDiv> *ref_to_geo_divs();
  boost::multi_array<int, 2> *ref_to_graticule_diagonals();
  boost::multi_array<XYPoint, 2> *ref_to_proj();
  FTReal2d *ref_to_rho_ft();
  FTReal2d *ref_to_rho_init();
  std::vector<std::vector<intersection> >
    scanlines_parallel_to_axis(bool, unsigned int) const;
  void set_area_errors();
  void set_geo_divs(const std::vector<GeoDiv>);
  void set_grid_dimensions(const unsigned int lx, const unsigned int ly);
  void set_inset_name(const std::string);
  void set_map_scale(const double);
  void set_pos(const std::string);
  void set_xmin(const unsigned int);
  void set_ymin(const unsigned int);
  bool target_area_is_missing(const std::string) const;
  double target_areas_at(const std::string) const;
  void target_areas_insert(const std::string, const double);
  void target_areas_replace(const std::string, const double);
  double total_inset_area() const;
  double total_target_area() const;

  // Write all intersections found to an EPS file named "*_intersections_*.eps"
  void write_intersections_to_eps(unsigned int);
};

#endif
