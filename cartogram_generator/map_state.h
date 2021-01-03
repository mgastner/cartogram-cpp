#ifndef MAP_STATE_H_
#define MAP_STATE_H_

#include "ft_real_2d.h"
#include "geo_div.h"
#include "colors.h"
#include <fftw3.h>
#include <vector>
#include <boost/multi_array.hpp>

struct XYPoint{ // For use in proj & proj2 arrays and flatten_density
  double x;
  double y;
};

struct XYPoint {
  double x;
  double y;
};

class MapState {
private:
  std::vector<GeoDiv> geo_divs_;
  std::map<std::string, double> target_areas;
  std::map<std::string, Color> colors;
  std::string id_header_;
  std::string visual_variable_file_;
  std::set<std::string> ids_in_visual_variables_file_;
  bool is_world_map_;
  bool write_density_to_eps_;
  unsigned int lx_, ly_;  // Lattice dimensions
  FTReal2d rho_init_;  // Rasterized density
  FTReal2d rho_ft_;  // Fourier transform
  fftw_plan fwd_plan_for_rho_, bwd_plan_for_rho_;
  unsigned int n_finished_integrations_;
  boost::multi_array<XYPoint, 2> proj_;
  boost::multi_array<XYPoint, 2> proj2_;
  MapState();
public:
  explicit MapState(const std::string, const bool, const bool);
  ~MapState();
  unsigned int n_geo_divs() const;
  const std::vector<GeoDiv> geo_divs() const;
  std::vector<GeoDiv> *ref_to_geo_divs();
  void target_areas_insert(std::string, double);
  void colors_insert(std::string, std::string);
  double target_areas_at(const std::string);
  const Color colors_at(const std::string);
  bool colors_empty() const;
  void set_id_header(const std::string);
  const std::string id_header() const;
  const std::string visual_variable_file() const;
  void insert_id_in_visual_variables_file(const std::string);
  const std::set<std::string> ids_in_visual_variables_file() const;
  bool is_world_map() const;
  bool trigger_write_density_to_eps() const;
  void make_grid(const unsigned int, const unsigned int);
  unsigned int lx() const;
  unsigned int ly() const;
  FTReal2d *ref_to_rho_init();
  FTReal2d *ref_to_rho_ft();
  void execute_fwd_plan() const;
  void execute_bwd_plan() const;
  void push_back(const GeoDiv);
  unsigned int n_finished_integrations() const;
  boost::multi_array<XYPoint, 2> *proj();
  boost::multi_array<XYPoint, 2> *proj2();
};

#endif
