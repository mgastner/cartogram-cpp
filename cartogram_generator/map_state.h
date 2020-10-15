#ifndef MAP_STATE_H_
#define MAP_STATE_H_

#include "geo_div.h"
#include <vector>

class MapState {
private:
  bool world;
  std::vector<GeoDiv> geo_divs;
  int lx, ly;  // Lattice dimensions
public:
  MapState(const bool);
  int n_geo_divs(void) const;
  std::vector<GeoDiv> get_geo_divs(void) const;
  std::vector<GeoDiv> *ref_to_geo_divs(void);
  bool is_world_map(void) const;
  void set_lx(const int);
  void set_ly(const int);
  void push_back(const GeoDiv);
};

#endif
