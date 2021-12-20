#include "cartogram_info.h"

CartogramInfo::CartogramInfo(const bool w,
                             const std::string v) :
  is_world_map_(w),
  visual_variable_file_(v)
{
  return;
}

double CartogramInfo::cart_total_target_area() const
{
  double area = 0.0;

  // Loop over inset states. Syntax from:
  // https://stackoverflow.com/questions/13087028/can-i-easily-iterate-over-
  // the-values-of-a-map-using-a-range-based-for-loop
  for (const auto &inset_info : inset_states_) {
       auto &inset_state = inset_info.second;
    area += inset_state.total_target_area();
  }
  return area;
}

void CartogramInfo::gd_to_inset_insert(const std::string id,
                                       const std::string inset)
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

const std::string CartogramInfo::inset_at_gd(const std::string id) const
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

bool CartogramInfo::original_ext_ring_is_clockwise() const
{
  return original_ext_ring_is_clockwise_;
}

std::map<std::string, InsetState> *CartogramInfo::ref_to_inset_states()
{
  return &inset_states_;
}

void CartogramInfo::replace_missing_and_zero_target_areas()
{
  // Check whether target areas exist that are missing or equal to zero
  bool ta_zero_exists = false;
  bool ta_na_exists = false;
  for (auto &inset_info : inset_states_) {
       auto &inset_state = inset_info.second;
    for (const auto &gd : inset_state.geo_divs()) {
      const double target_area = inset_state.target_areas_at(gd.id());
      if (target_area < 0.0) {
        ta_na_exists = true;

        // Insert true into an std::unordered_map that stores whether a
        // GeoDiv's target area is missing
        inset_state.is_input_target_area_missing_insert(gd.id(), true);
      } else if (target_area == 0) {
        ta_zero_exists = true;
      }

      // Insert false if the value does not already exist in the map
      inset_state.is_input_target_area_missing_insert(gd.id(), false);
    }
  }

  // Deal with target areas that are equal to zero
  if (ta_zero_exists) {

    // Find smallest positive target area
    double min_positive_area = dbl_inf;
    for (const auto &inset_info : inset_states_) {
      auto &inset_state = inset_info.second;
      for (const auto &gd : inset_state.geo_divs()) {
        const double target_area = inset_state.target_areas_at(gd.id());

        // Target area not equal to zero and not missing
        if (target_area > 0.0) {
          min_positive_area = std::min(min_positive_area, target_area);
        }
      }
    }

    // If all target areas are zero or missing, we assign the minimum GeoDiv
    // area (instead of the minimum target area) to min_positive_area
    if (min_positive_area == dbl_inf) {
      for (const auto &inset_info : inset_states_) {
           auto &inset_state = inset_info.second;
        for (auto const &gd : inset_state.geo_divs()) {
          min_positive_area = std::min(min_positive_area, gd.area());
        }
      }
    }

    // Replace non-positive target areas with a fraction of the smallest
    // positive target area
    const double replacement_for_nonpositive_area = 0.1 * min_positive_area;
    std::cerr << "Replacing zero target area with "
              << replacement_for_nonpositive_area
              << " (0.1 times the minimum positive area)."
              << std::endl;
    for (auto &inset_info : inset_states_) {
      auto &inset_state = inset_info.second;
      for (const auto &gd : inset_state.geo_divs()) {
        const double target_area = inset_state.target_areas_at(gd.id());
        if (target_area == 0.0) {
          inset_state.target_areas_replace(gd.id(),
                                           replacement_for_nonpositive_area);
        }
      }
    }
  }

  // Deal with missing target areas
  if (ta_na_exists) {
    double total_cart_non_na_area = 0.0;
    double total_cart_non_na_ta = 0.0;

    for (const auto &inset_info : inset_states_) {
      auto &inset_state = inset_info.second;
      for (const auto &gd : inset_state.geo_divs()) {
        if (!inset_state.target_area_is_missing(gd.id())) {
          total_cart_non_na_area += gd.area();
          total_cart_non_na_ta += inset_state.target_areas_at(gd.id());
        }
      }
    }

    // Assign new target areas to GeoDivs
    for (auto &inset_info : inset_states_) {
      auto &inset_state = inset_info.second;
      for (const auto &gd : inset_state.geo_divs()) {
        if (inset_state.target_area_is_missing(gd.id())) {
          double new_target_area;

          // If all target areas are missing, make all GeoDivs equal to their
          // geographic area
          if (total_cart_non_na_ta == 0.0) {
            new_target_area = gd.area();
          } else {

            // Replace target_area
            new_target_area =
              (total_cart_non_na_ta / total_cart_non_na_area) * gd.area();
          }
          inset_state.target_areas_replace(gd.id(), new_target_area);
        }
      }
    }
  }
}

void CartogramInfo::set_id_header(const std::string id)
{
  id_header_ = id;
  return;
}

void CartogramInfo::set_original_ext_ring_is_clockwise(
  const bool original_ext_ring_is_clockwise)
{
  original_ext_ring_is_clockwise_ = original_ext_ring_is_clockwise;
  return;
}

const std::string CartogramInfo::visual_variable_file() const
{
  return visual_variable_file_;
}
