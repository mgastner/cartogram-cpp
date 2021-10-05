#ifndef INSET_STATE_H_
#define INSET_STATE_H_

#include "ft_real_2d.h"
#include "geo_div.h"
#include "colors.h"
#include <vector>
#include <boost/multi_array.hpp>
#include <map>

// Struct to store the X and Y coordinates of a 2D point
struct XYPoint {
  double x;
  double y;

  // Overload "<" operator for this data type. Idea from
  // https://stackoverflow.com/questions/4892680/sorting-a-vector-of-structs
  // To be able to sort XYPoint
  bool operator < (const XYPoint &rhs) const
  {
    if (x != rhs.x) {
      return (x < rhs.x);
    }
    return y < rhs.y;
  }

  // Overload "==" and "!=" operators for XYPoint type.
  bool operator == (const XYPoint &rhs) const
  {
    return (x == rhs.x && y == rhs.y);
  }
};

// Struct to store intersection data
struct intersection {
  double x;  // x-coordinate of intersection
  double target_density;  // GeoDiv's target_density
  std::string geo_div_id;  // GeoDIv's ID
  bool direction;  // Does intersection enter (true) or exit (false)?

  // Overloading "<" operator, similar to above
  bool operator < (const intersection &rhs) const
  {
    return (x < rhs.x || (x == rhs.x && direction < rhs.direction));
  }
};

struct max_area_error_info {
  double value;
  std::string geo_div;
};

class InsetState {
private:
  std::unordered_map<std::string, double> area_errors_;
  CGAL::Bbox_2 bbox_;  // Bounding box
  std::unordered_map<std::string, Color> colors_;
  fftw_plan fwd_plan_for_rho_, bwd_plan_for_rho_;
  std::vector<GeoDiv> geo_divs_;  // Geographic divisions in this inset

  // Chosen diagonal for each graticule cell
  boost::multi_array<int, 2> graticule_diagonals_;
  // Horizontal and vertical adjacency graphs
  std::vector<std::vector<intersection> > horizontal_adj_, vertical_adj_;
  std::string inset_name_; // Map name, appended with Position if n_insets > 2
  unsigned int lx_, ly_;  // Lattice dimensions
  double map_scale_; // Double to map scale
  unsigned int new_xmin_, new_ymin_; // To store map translation vector
  unsigned int n_finished_integrations_;
  std::string pos_;  // Position of inset ("C", "T" etc.)
  boost::multi_array<XYPoint, 2> proj_;  // Cartogram projection

  // Rasterized density and its Fourier transform
  FTReal2d rho_init_, rho_ft_;
  std::unordered_map<std::string, double> target_areas_;

  // Make default contructor private so that only
  // InsetState(const std::string) can be called as constructor
  InsetState();
public:
  explicit InsetState(const std::string);  // Constructor
  double area_errors_at(const std::string) const;
  CGAL::Bbox_2 bbox() const;
  double cart_area() const;
  bool color_found(const std::string id) const;
  const Color colors_at(const std::string) const;
  bool colors_empty() const;
  void colors_insert(const std::string, const Color);
  void colors_insert(const std::string, std::string);
  unsigned int colors_size() const;
  void destroy_fftw_plans_for_rho();
  void execute_fftw_bwd_plan() const;
  void execute_fftw_fwd_plan() const;
  const std::vector<GeoDiv> geo_divs() const;
  boost::multi_array<int, 2> *graticule_diagonals();
  const std::vector<std::vector<intersection> > horizontal_adj() const;
  void increment_integration();
  const std::string inset_name() const;
  unsigned int lx() const;
  unsigned int ly() const;
  void make_fftw_plans_for_rho();
  double map_scale() const;
  double max_area_error_stdcerr() const;
  struct max_area_error_info max_area_error() const;
  unsigned int new_xmin() const;
  unsigned int new_ymin() const;
  unsigned int n_finished_integrations() const;
  unsigned int n_geo_divs() const;
  const std::string pos() const;
  boost::multi_array<XYPoint, 2> *proj();
  void push_back(const GeoDiv);
  std::vector<GeoDiv> *ref_to_geo_divs();
  FTReal2d *ref_to_rho_ft();
  FTReal2d *ref_to_rho_init();
  void set_area_errors();
  void set_geo_divs(std::vector<GeoDiv>);
  void set_grid_dimensions(unsigned int lx, unsigned int ly);
  void set_horizontal_adj(std::vector<std::vector<intersection> >);
  void set_inset_name(std::string);
  void set_map_scale(const double);
  void set_pos(std::string);
  void set_vertical_adj(std::vector<std::vector<intersection> >);
  void set_xmin(const unsigned int);
  void set_ymin(const unsigned int);
  bool target_area_is_missing(const std::string) const;
  double target_areas_at(const std::string) const;
  void target_areas_insert(std::string, double);
  void target_areas_replace(std::string, double);
  double total_target_area() const;
  const std::vector<std::vector<intersection> > vertical_adj() const;
};

#endif
