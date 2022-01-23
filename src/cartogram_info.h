#ifndef CARTOGRAM_INFO_H_
#define CARTOGRAM_INFO_H_

#include "inset_state.h"
#include "constants.h"
#include <vector>

class CartogramInfo {
private:
  std::map<std::string, std::string> gd_to_inset_;
  std::string id_header_;
  std::set<std::string> ids_in_visual_variables_file_;
  std::map<std::string, InsetState> inset_states_;

  // TODO: We assume that either all external rings are counterclockwise or
  // all are clockwise. This dichotomy covers most geospatial boundary files
  // in the wild, but it would still be sensible to allow cases where there
  // are external rings with opposite winding directions.
  bool original_ext_ring_is_clockwise_;
  bool is_world_map_;
  std::string visual_variable_file_;

public:
  explicit CartogramInfo(const bool, const std::string);
  double cart_total_target_area() const;
  void gd_to_inset_insert(const std::string, const std::string);
  const std::string id_header() const;
  const std::set<std::string> ids_in_visual_variables_file() const;
  void insert_id_in_visual_variables_file(const std::string);
  void insert_inset_state(const std::string, const InsetState);
  const std::string inset_at_gd(const std::string) const;
  const std::map<std::string, InsetState> inset_states() const;
  bool is_world_map() const;
  unsigned int n_insets() const;
  bool original_ext_ring_is_clockwise() const;
  std::map<std::string, InsetState> *ref_to_inset_states();
  void replace_missing_and_zero_target_areas();
  void set_id_header(const std::string);
  void set_original_ext_ring_is_clockwise(const bool);
  const std::string visual_variable_file() const;
};
#endif
