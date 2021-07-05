#include "cartogram_info.h"

CartogramInfo::CartogramInfo(std::string v, const bool w, const bool wd2eps) :
  visual_variable_file_(v),
  is_world_map_(w),
  write_density_to_eps_(wd2eps)
{
  return;
}

void CartogramInfo::gd_to_inset_insert(const std::string id,
                                       std::string inset)
{
  gd_to_inset_.insert(std::pair<std::string, std::string>(id, inset));
  return;
}

const std::string CartogramInfo::id_header() const
{
  return id_header_;
}

const std::set<std::string> CartogramInfo::ids_in_visual_variables_file()
const
{
  return ids_in_visual_variables_file_;
}

void CartogramInfo::insert_id_in_visual_variables_file(const std::string id)
{
  ids_in_visual_variables_file_.insert(id);
}

const std::string CartogramInfo::inset_at_gd(const std::string id)
{
  return gd_to_inset_.at(id);
}

const std::vector<InsetState> CartogramInfo::inset_states() const
{
  return inset_states_;
}

bool CartogramInfo::is_world_map() const
{
  return is_world_map_;
}

unsigned int CartogramInfo::n_insets() const
{
  return inset_states_.size();
}

void CartogramInfo::push_back(const InsetState is)
{
  inset_states_.push_back(is);
  return;
}

std::vector<InsetState> *CartogramInfo::ref_to_inset_states()
{
  return &inset_states_;
}

void CartogramInfo::set_id_header(const std::string id)
{
  id_header_ = id;
  return;
}

void CartogramInfo::set_total_cart_target_area(double total_cart_target_area)
{
  total_cart_target_area_ = total_cart_target_area;
}

double CartogramInfo::total_cart_target_area()
{
  return total_cart_target_area_;
}

bool CartogramInfo::trigger_write_density_to_eps() const
{
  return write_density_to_eps_;
}

const std::string CartogramInfo::visual_variable_file() const
{
  return visual_variable_file_;
}
