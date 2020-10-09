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
  GeoDiv get_geo_div(const unsigned int);
  bool is_world_map(void);
};

#endif
