#ifndef MAP_STATE_H_
#define MAP_STATE_H_

#include "ft_real_2d.h"
#include "geo_div.h"
#include <fftw3.h>
#include <vector>

class MapState {
private:
  bool world;
  std::vector<GeoDiv> geo_divs;
  unsigned int lx = 0, ly = 0;  // Lattice dimensions
  FTReal2d rho_init;  // Rasterized density
  FTReal2d rho_ft;  // Fourier transform
  fftw_plan plan_fwd, plan_bwd;  // Plan the Fourier transform
  MapState();
public:
  explicit MapState(const bool);
  ~MapState();
  const unsigned int n_geo_divs() const;
  const std::vector<GeoDiv> get_geo_divs() const;
  std::vector<GeoDiv> *ref_to_geo_divs();
  const bool is_world_map() const;
  void make_grid(const unsigned int, const unsigned int);
  const unsigned int get_lx() const;
  const unsigned int get_ly() const;
  FTReal2d *ref_to_rho_init();
  FTReal2d *ref_to_rho_ft();
  const fftw_plan get_plan_fwd() const;
  const fftw_plan get_plan_bwd() const;
  void free_rho();
  void push_back(const GeoDiv);
};

#endif
