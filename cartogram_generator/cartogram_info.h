#ifndef CARTOGRAM_INFO_H_
#define CARTOGRAM_INFO_H_

#include "inset_state.h"
#include <vector>

class CartogramInfo {
private:
  std::vector<InsetState> inset_states_;
  std::map<std::string, std::string> gd_to_inset_;
  std::string id_header_;
  std::string visual_variable_file_;
  std::set<std::string> ids_in_visual_variables_file_;
  bool is_world_map_;
  bool write_density_to_eps_;
  std::string map_name_;
  double total_cart_target_area_;

  // To store all bbox with position
  std::map <std::string, CGAL::Bbox_2> all_bbox_with_pos_; 

  // To store all frame bbox with position
  std::map <std::string, CGAL::Bbox_2> all_frame_bbox_with_pos_; 


public:
  explicit CartogramInfo(const std::string, const bool, const bool);
  void set_id_header(const std::string);
  void insert_id_in_visual_variables_file(const std::string);
  const std::set<std::string> ids_in_visual_variables_file() const;
  const std::string id_header() const;
  const std::string visual_variable_file() const;
  bool is_world_map() const;
  bool trigger_write_density_to_eps() const;
  void set_map_name(std::string map_name);
  unsigned int n_insets() const;
  const std::vector<InsetState> inset_states() const;
  std::vector<InsetState> *ref_to_inset_states();
  void push_back(const InsetState);
  void gd_to_inset_insert(std::string, std::string);
  const std::string inset_at_gd(const std::string);

  // Stores total target area of all GeoDivs
  void set_total_cart_target_area(double); 

  // Returns total target area of all GeoDivs
  double total_cart_target_area();

  // Stores all bbox with position
  void set_bbox_at_pos(std::string, CGAL::Bbox_2); 

  // Returns all bbox with position
  std::map <std::string, CGAL::Bbox_2> all_bbox_with_pos(); 

  // Stores frame bbox with position
  void set_frame_bbox_at_pos(std::string, CGAL::Bbox_2); 

  // Returns all frame bbox with position
  std::map <std::string, CGAL::Bbox_2> all_frame_bbox_with_pos(); 
};
#endif
