#ifndef MAP_STATE_H_
#define MAP_STATE_H_

#include "geo_div.h"
#include <fftw3.h>
#include <vector>

class MapState {
private:
  bool world;
  std::vector<GeoDiv> geo_divs;
  int lx, ly;  // Lattice dimensions
  double *rho_init = NULL;  // Raterized density
  double *rho_ft = NULL;  // Fourier transform
  fftw_plan plan_fwd, plan_bwd;  // Plan the Fourier transform
  MapState();
public:
  explicit MapState(const bool);
  ~MapState();
  int n_geo_divs(void) const;
  std::vector<GeoDiv> get_geo_divs(void) const;
  std::vector<GeoDiv> *ref_to_geo_divs(void);
  bool is_world_map(void) const;
  void make_grid(const unsigned int, const unsigned int);
  int get_lx(void);
  int get_ly(void);
  double *get_rho_init(void);
  double *get_rho_ft(void);
  fftw_plan get_plan_fwd(void);
  void free_rho(void);
  void push_back(const GeoDiv);
};

#endif
