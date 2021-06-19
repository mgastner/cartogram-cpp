#ifndef INSET_STATE_H_
#define INSET_STATE_H_

#include "ft_real_2d.h"
#include "geo_div.h"
#include "colors.h"
#include <fftw3.h>
#include <vector>
#include <boost/multi_array.hpp>

// Struct to store the X and Y coordinates of a 2D point
struct XYPoint{
  double x;
  double y;
};

// Struct to store intersection data
struct intersection {
  double x;  // x-coordinate of intersection
  double target_density;  // GeoDiv's target_density
  std::string geo_div_id; // GeoDIv's ID
  bool direction;  // Does intersection enter (true) or exit (false)?

  // Overload "<" operator for this data type. Idea from
  // https://stackoverflow.com/questions/4892680/sorting-a-vector-of-structs
  bool operator < (const intersection &rhs) const
  {
    return (x < rhs.x || (x == rhs.x && direction < rhs.direction));
  }
};

class InsetState {
private:
  std::string pos_; // Position of inset ("C", "TR" etc.)
  std::string inset_name_; // Map name, appended with Position if n_insets > 2
  std::vector<GeoDiv> geo_divs_;
  std::unordered_map<std::string, double> target_areas;
  std::unordered_map<std::string, double> area_errs;
  std::unordered_map<std::string, Color> colors;
  unsigned int lx_, ly_;  // Lattice dimensions
  unsigned int new_xmin_, new_ymin_; // To store map translation vector
  double map_scale_; // Double to map scale
  FTReal2d rho_init_;  // Rasterized density
  FTReal2d rho_ft_;  // Fourier transform
  fftw_plan fwd_plan_for_rho_, bwd_plan_for_rho_;
  unsigned int n_finished_integrations_;
  boost::multi_array<XYPoint, 2> proj_;

  // Horizontal and Vertical adjacency graphs
  std::vector<std::vector<intersection> > horizontal_adj_;
  std::vector<std::vector<intersection> > vertical_adj_;
  InsetState();
public:
  explicit InsetState(const std::string);
  ~InsetState();
  unsigned int n_geo_divs() const;
  const std::vector<GeoDiv> geo_divs() const;
  std::vector<GeoDiv> *ref_to_geo_divs();
  void set_geo_divs(std::vector<GeoDiv>);
  void target_areas_insert(std::string, double);
  void colors_insert(const std::string, std::string);
  void colors_insert(const std::string, const Color);
  double target_areas_at(const std::string);
  bool target_area_is_missing(const std::string) const;
  const Color colors_at(const std::string);
  bool colors_empty() const;
  bool color_found(const std::string id) const;
  void make_grid(const unsigned int, const unsigned int);
  unsigned int lx() const;
  unsigned int ly() const;
  unsigned int new_xmin() const;
  unsigned int new_ymin() const;
  void set_new_xmin(const unsigned int);
  void set_new_ymin(const unsigned int);
  double map_scale() const;
  void set_map_scale(const double);
  FTReal2d *ref_to_rho_init();
  FTReal2d *ref_to_rho_ft();
  void execute_fwd_plan() const;
  void execute_bwd_plan() const;
  void push_back(const GeoDiv);
  unsigned int n_finished_integrations() const;
  void inc_integration();
  boost::multi_array<XYPoint, 2> *proj();
  void set_area_errs();
  double area_errs_at(const std::string) const;
  double max_area_err() const;
  void set_pos(std::string);
  const std::string pos() const;
  void set_inset_name(std::string);
  const std::string inset_name() const;
  void set_horizontal_adj(std::vector<std::vector<intersection> >);
  void set_vertical_adj(std::vector<std::vector<intersection> >);
  const std::vector<std::vector<intersection> > horizontal_adj() const;
  const std::vector<std::vector<intersection> > vertical_adj() const;
  unsigned int colors_size() const;
};

#endif
