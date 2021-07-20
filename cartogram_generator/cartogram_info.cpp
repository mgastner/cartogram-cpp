#include "cartogram_info.h"

CartogramInfo::CartogramInfo(const bool w,
                             const std::string v,
                             const bool wd2eps) :
  is_world_map_(w),
  visual_variable_file_(v),
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

void CartogramInfo::insert_inset_state(const std::string inset_pos,
                                       const InsetState inset_state)
{
  inset_states_.insert(std::pair<std::string, InsetState>(inset_pos,
                                                          inset_state));
  return;
}

const std::string CartogramInfo::inset_at_gd(const std::string id)
{
  return gd_to_inset_.at(id);
}

const std::map<std::string, InsetState> CartogramInfo::inset_states() const
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

std::map<std::string, InsetState> *CartogramInfo::ref_to_inset_states()
{

  return &inset_states_;
}

void CartogramInfo::set_id_header(const std::string id)
{
  id_header_ = id;
  return;
}

double CartogramInfo::total_cart_target_area() const
{
  double area = 0.0;

  // Loop over inset states. Syntax from:
  // https://stackoverflow.com/questions/13087028/can-i-easily-iterate-over-
  // the-values-of-a-map-using-a-range-based-for-loop
  for (auto const &inset_state : inset_states_ | std::views::values) {
    for (auto gd : inset_state.geo_divs()) {
      area += inset_state.target_areas_at(gd.id());
    }
  }
  return area;
}

bool CartogramInfo::trigger_write_density_to_eps() const
{
  return write_density_to_eps_;
}

const std::string CartogramInfo::visual_variable_file() const
{
  return visual_variable_file_;
}
