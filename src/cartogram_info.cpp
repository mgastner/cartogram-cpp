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

const std::string CartogramInfo::id_header() const
{
  return id_header_;
}

const std::set<std::string> CartogramInfo::ids_in_visual_variables_file()
const
{
  return ids_in_visual_variables_file_;
}

void CartogramInfo::insert_gd_into_inset(const std::string id,
                                         const std::string inset)
{
  gd_to_inset_.insert(std::pair<std::string, std::string>(id, inset));
  return;
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

  // Get total current area and total target area
  double total_start_area_with_data = 0.0;
  double total_target_area_with_data = 0.0;
  for (const auto &inset_info : inset_states_) {
    auto &inset_state = inset_info.second;
    for (const auto &gd : inset_state.geo_divs()) {
      if (!inset_state.target_area_is_missing(gd.id())) {
        total_start_area_with_data += gd.area();
        total_target_area_with_data += inset_state.target_area_at(gd.id());
      }
    }
  }

  // Calculate threshold for small target areas. For GeoDivs below this
  // threshold, the target area is scaled up for easier calculation.
  const double small_target_area_threshold =
    total_target_area_with_data * small_area_threshold_frac;

  // Check whether target areas exist that are missing or very small
  bool small_target_area_exists = false;
  bool missing_target_area_exists = false;
  for (auto &inset_info : inset_states_) {
    auto &inset_state = inset_info.second;
    for (const auto &gd : inset_state.geo_divs()) {
      const double target_area = inset_state.target_area_at(gd.id());
      inset_state.insert_whether_input_target_area_is_missing(
        gd.id(),
        target_area < 0.0);
      if (target_area < 0.0) {
        missing_target_area_exists = true;
      } else if (target_area <= small_target_area_threshold) {
        small_target_area_exists = true;
      }
    }
  }

  // Deal with target areas that are below the threshold
  if (small_target_area_exists) {
    double replacement_target_area;

    // We replace the zero and small areas, if any, with
    // small_area_absolute_threshold if not all target areas are initially
    // missing or zero
    if (small_target_area_threshold > 0.0) {
      std::cerr << "Replacing small target areas."
                << std::endl;
      replacement_target_area = small_target_area_threshold;
    } else {

      // If all target areas are zero or missing, we assign the minimum GeoDiv
      // area (instead of the minimum target area) as replacement_target_area.
      std::cerr << "No non-zero target area.\n"
                << "Setting zero target areas to the minimum positive area."
                << std::endl;
      double min_positive_area = dbl_inf;
      for (const auto &inset_info : inset_states_) {
        auto &inset_state = inset_info.second;
        for (const auto &gd : inset_state.geo_divs()) {
          min_positive_area = std::min(min_positive_area, gd.area());
        }
      }
      replacement_target_area = min_positive_area;
    }

    // Replace the small target areas
    for (auto &inset_info : inset_states_) {
      auto &inset_state = inset_info.second;
      for (const auto &gd : inset_state.geo_divs()) {

        // Current target area
        const double target_area = inset_state.target_area_at(gd.id());
        if ((target_area >= 0.0) &&
            (target_area <= small_target_area_threshold)) {
          inset_state.replace_target_area(gd.id(), replacement_target_area);
          std::cerr << gd.id() << ": "
                    << target_area << " to " << replacement_target_area
                    << std::endl;

          // Update total target area
          total_target_area_with_data +=
            (replacement_target_area - target_area);
        }
      }
    }
  }

  // Deal with missing target areas
  if (missing_target_area_exists) {

    // Assign new target areas to GeoDivs
    for (auto &inset_info : inset_states_) {
      auto &inset_state = inset_info.second;
      for (const auto &gd : inset_state.geo_divs()) {
        if (inset_state.target_area_is_missing(gd.id())) {
          double new_target_area;

          // If all target areas are missing, make all GeoDivs equal to their
          // geographic area
          if (total_target_area_with_data == 0.0) {
            new_target_area = gd.area();
          } else {

            // Replace target_area
            const double mean_density =
              (total_target_area_with_data / total_start_area_with_data);
            new_target_area = mean_density * gd.area();
          }
          inset_state.replace_target_area(gd.id(), new_target_area);
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

void CartogramInfo::set_map_name(const std::string map_name)
{
  map_name_ = map_name;
  return;
}

const std::string CartogramInfo::map_name() const
{
  return map_name_;
}
