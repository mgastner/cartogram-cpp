#ifndef CARTOGRAM_INFO_H_
#define CARTOGRAM_INFO_H_

#include "inset_state.h"
#include <vector>

class CartogramInfo {
private:
  std::map<std::string, std::string> gd_to_inset_;
  std::string id_header_;
  std::set<std::string> ids_in_visual_variables_file_;
  std::vector<InsetState> inset_states_;
  bool is_world_map_;
  double total_cart_target_area_;
  std::string visual_variable_file_;
  bool write_density_to_eps_;
public:
  explicit CartogramInfo(const std::string, const bool, const bool);
  void gd_to_inset_insert(std::string, std::string);
  const std::string id_header() const;
  const std::set<std::string> ids_in_visual_variables_file() const;
  void insert_id_in_visual_variables_file(const std::string);
  const std::string inset_at_gd(const std::string);
  const std::vector<InsetState> inset_states() const;
  bool is_world_map() const;
  unsigned int n_insets() const;
  void push_back(const InsetState);
  std::vector<InsetState> *ref_to_inset_states();
  void set_id_header(const std::string);
  void set_total_cart_target_area(double);
  double total_cart_target_area();
  bool trigger_write_density_to_eps() const;
  const std::string visual_variable_file() const;
};
#endif
