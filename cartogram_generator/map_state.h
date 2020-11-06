#ifndef MAP_STATE_H_
#define MAP_STATE_H_

#include "ft_real_2d.h"
#include "geo_div.h"
#include <fftw3.h>
#include <vector>

class MapState {
private:
  bool is_world_map_;
  std::vector<GeoDiv> geo_divs_;
  std::map<std::string, double> target_areas;
  std::map<std::string, std::string> colors;
  unsigned int lx_ = 0, ly_ = 0;  // Lattice dimensions
  FTReal2d rho_init_;  // Rasterized density
  FTReal2d rho_ft_;  // Fourier transform
  fftw_plan fwd_plan_, bwd_plan_;  // Plan the Fourier transform
  MapState();
public:
  explicit MapState(const bool);
  ~MapState();
  const unsigned int n_geo_divs() const;
  const std::vector<GeoDiv> geo_divs() const;
  std::vector<GeoDiv> *ref_to_geo_divs();
  void target_areas_insert(std::string, double);
  void colors_insert(std::string, std::string);
  const double target_areas_at(const std::string id);
  const std::string colors_at(const std::string id);
  const bool is_world_map() const;
  void make_grid(const unsigned int, const unsigned int);
  const unsigned int lx() const;
  const unsigned int ly() const;
  FTReal2d *ref_to_rho_init();
  FTReal2d *ref_to_rho_ft();
  void execute_fwd_plan() const;
  void execute_bwd_plan() const;
  void push_back(const GeoDiv);
};

#endif
