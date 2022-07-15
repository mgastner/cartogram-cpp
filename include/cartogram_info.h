#ifndef CARTOGRAM_INFO_H_
#define CARTOGRAM_INFO_H_

#include "inset_state.h"
#include "argparse.hpp"
#include "constants.h"
#include <vector>

class CartogramInfo {
private:
  std::map<std::string, std::string> gd_to_inset_;
  std::string id_header_;
  std::set<std::string> ids_in_visual_variables_file_;
  std::map<std::string, InsetState> inset_states_;
  bool is_world_map_;
  std::string map_name_;

  // TODO: We assume that either all external rings are counterclockwise or
  // all are clockwise. This dichotomy covers most geospatial boundary files
  // in the wild, but it would still be sensible to allow cases where there
  // are external rings with opposite winding directions.
  bool original_ext_ring_is_clockwise_;
  std::string visual_variable_file_;
  nlohmann::json cgal_to_json();

public:
  explicit CartogramInfo(const bool, const std::string);
  double cart_total_target_area() const;
  double area() const;
  bool is_world_map() const;
  const std::string map_name() const;
  unsigned int n_geo_divs() const;
  unsigned int n_insets() const;
  bool original_ext_ring_is_clockwise() const;
  void read_csv(argparse::ArgumentParser);
  void read_geojson(const std::string, const bool, std::string*);
  std::map<std::string, InsetState> *ref_to_inset_states();
  void replace_missing_and_zero_target_areas();
  void set_map_name(const std::string);
  void shift_insets_to_target_position();
  void write_geojson(std::string, std::string, std::ostream &, bool);
};
#endif
