#ifndef INSET_STATE_H_
#define INSET_STATE_H_

#include "ft_real_2d.h"
#include "geo_div.h"
#include "colors.h"
#include <fftw3.h>
#include <vector>
#include <boost/multi_array.hpp>

struct XYPoint{
  double x;
  double y;
};

class InsetState {
private:
  std::string inset_pos_;
  std::vector<GeoDiv> geo_divs_;
  std::map<std::string, double> target_areas;
  std::map<std::string, double> area_errs;
  std::map<std::string, Color> colors;
  unsigned int lx_, ly_;  // Lattice dimensions
  unsigned int new_xmin_, new_ymin_; // To store map translation vector
  double map_scale_; // Double to map scale
  FTReal2d rho_init_;  // Rasterized density
  FTReal2d rho_ft_;  // Fourier transform
  fftw_plan fwd_plan_for_rho_, bwd_plan_for_rho_;
  unsigned int n_finished_integrations_;
  boost::multi_array<XYPoint, 2> proj_;
  InsetState();
public:
  explicit InsetState(const std::string);
  ~InsetState();
  unsigned int n_geo_divs() const;
  const std::vector<GeoDiv> geo_divs() const;
  std::vector<GeoDiv> *ref_to_geo_divs();
  void set_geo_divs(std::vector<GeoDiv>);
  void target_areas_insert(std::string, double);
  void colors_insert(std::string, std::string);
  double target_areas_at(const std::string);
  bool target_area_is_missing(const std::string) const;
  const Color colors_at(const std::string);
  bool colors_empty() const;
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
  void set_inset_pos(std::string inset_pos);
  const std::string inset_pos() const;
};

#endif
