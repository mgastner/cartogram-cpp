#include "cartogram_info.h"

CartogramInfo::CartogramInfo(const bool w,
                             const std::string v) :
  is_world_map_(w),
  visual_variable_file_(v)
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

bool CartogramInfo::original_ext_ring_is_clockwise()
{
  return original_ext_ring_is_clockwise_;
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

double CartogramInfo::cart_non_missing_target_area() const
{
  double area = 0.0;

  // Loop over inset states. Syntax from:
  // https://stackoverflow.com/questions/13087028/can-i-easily-iterate-over-
  // the-values-of-a-map-using-a-range-based-for-loop
  for (auto const &inset_state : inset_states_ | std::views::values) {
    for (auto gd : inset_state.geo_divs()) {
      if (!inset_state.target_area_is_missing(gd.id())) {
        area += inset_state.target_areas_at(gd.id());
      }
    }
  }
  return area;
}

const std::string CartogramInfo::visual_variable_file() const
{
  return visual_variable_file_;
}

void CartogramInfo::set_original_ext_ring_is_clockwise(
  bool original_ext_ring_is_clockwise)
{
  original_ext_ring_is_clockwise_ = original_ext_ring_is_clockwise;
  return;
}

void CartogramInfo::replace_missing_and_zero_target_areas()
{


  // Checking whether target areas that are equal to zero or NA exist
  bool ta_zero_exists = false, ta_na_exists = false;
  for (auto &inset_state : inset_states_ | std::views::values) {
    for (auto const gd : inset_state.geo_divs()) {
      double target_area = inset_state.target_areas_at(gd.id());
      if (target_area < 0.0) {
        ta_na_exists = true;

        // Inserting true into map which stores whether a GeoDivs has target
        // area missing
        inset_state.is_target_area_missing_insert(gd.id(), true);
      } else if (target_area == 0) {
        ta_zero_exists = true;
      }

      // Will only insert false if the value does not already exist in the map
      inset_state.is_target_area_missing_insert(gd.id(), false);
    }
  }

  // Dealing with target areas that are equal to zero
  if (ta_zero_exists) {

    // Find smallest positive target area
    double min_positive_area = dbl_inf;
    for (auto const &inset_state : inset_states_ | std::views::values) {
      for (auto const gd : inset_state.geo_divs()) {
        double target_area = inset_state.target_areas_at(gd.id());

        // Target area not equal to zero and not missing
        if (target_area > 0.0) {
          min_positive_area = std::min(min_positive_area, target_area);
        }
      }
    }

    // If all target areas are zero or missing, we assign the minimum GeoDiv
    // area (instead of the minimum target area) to min_positive_area
    if (min_positive_area == dbl_inf) {
      for (auto const &inset_state : inset_states_ | std::views::values) {
        for (auto const gd : inset_state.geo_divs()) {
          min_positive_area = std::min(min_positive_area, gd.area());
        }
      }
    }

    // Replace non-positive target areas with a fraction of the smallest
    // positive target area
    double replacement_for_nonpositive_area = 0.1 * min_positive_area;
    std::cerr << "Replacing zero target area with "
              << replacement_for_nonpositive_area
              << " (0.1 times the minimum positive area)."
              << std::endl;
    for (auto &inset_state : inset_states_ | std::views::values) {
      for (auto const gd : inset_state.geo_divs()) {
        double target_area = inset_state.target_areas_at(gd.id());
        if (target_area == 0.0) {
          inset_state.target_areas_replace(gd.id(),
                                           replacement_for_nonpositive_area);
        }
      }
    }
  }

  // Dealing with missing target areas
  if (ta_na_exists) {
    double total_non_na_area = 0.0;
    double total_non_na_ta = cart_non_missing_target_area();

    // Finding total area of non na geodivs
    for (auto const &inset_state : inset_states_ | std::views::values) {
      total_non_na_area += inset_state.cart_area();
    }

    // Assigning GeoDivs new target areas
    for (auto &inset_state : inset_states_ | std::views::values) {
      for (auto const gd : inset_state.geo_divs()) {
        if (inset_state.target_area_is_missing(gd.id())) {
          double new_target_area;

          // If all target areas are missing, make all GeoDivs equal to their
          // geographic area
          if (total_non_na_ta == 0.0) {
            new_target_area = gd.area();
          } else {

            // Replace target_area
            new_target_area = (total_non_na_ta / total_non_na_area) * gd.area();
          }

          inset_state.target_areas_replace(gd.id(), new_target_area);
        }
      }
    }
  }
}
