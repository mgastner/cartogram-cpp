#ifndef MAP_STATE_H_
#define MAP_STATE_H_

#include "geo_div.h"
#include <vector>

class MapState {
private:
  bool world;
  std::vector<GeoDiv> geo_divs;
public:
  MapState(const bool);
  void push_back(const GeoDiv);
  int n_geo_divs(void);
  std::vector<GeoDiv> get_geo_divs(void);
  bool is_world_map(void);
};

#endif
