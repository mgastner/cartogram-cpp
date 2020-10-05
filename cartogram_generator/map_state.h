#ifndef MAP_STATE_H_
#define MAP_STATE_H_

class MapState {
private:
  int n_geo_divs;
  bool world;
public:
  MapState(bool);
  int get_n_geo_divs(void);
  bool is_world_map(void);
};

#endif
