#include "cartogram_info.hpp"
#include "constants.hpp"
#include "csv.hpp"
#include <iostream>
#include <utility>

CartogramInfo::CartogramInfo(const bool w, const std::string &v)
    : is_world_map_(w), visual_variable_file_(std::move(v))
{
}

double CartogramInfo::cart_initial_total_target_area() const
{
  double target_area = 0.0;

  // Iterate over inset states. Syntax from:
  // https://stackoverflow.com/questions/13087028/can-i-easily-iterate-over-
  // the-values-of-a-map-using-a-range-based-for-loop
  for (const auto &inset_info : inset_states_) {
    auto &inset_state = inset_info.second;
    target_area += inset_state.initial_target_area();
  }
  return target_area;
}

double CartogramInfo::area() const
{
  double area = 0.0;

  // Iterate over inset states. Syntax from:
  // https://stackoverflow.com/questions/13087028/can-i-easily-iterate-over-
  // the-values-of-a-map-using-a-range-based-for-loop
  for (const auto &inset_info : inset_states_) {
    auto &inset_state = inset_info.second;
    area += inset_state.total_inset_area();
  }
  return area;
}

bool CartogramInfo::is_world_map() const
{
  return is_world_map_;
}

unsigned int CartogramInfo::n_geo_divs() const
{
  unsigned int n_geo_divs = 0;
  for (const auto &inset_info : inset_states_) {
    const auto &inset_state = inset_info.second;
    n_geo_divs += inset_state.n_geo_divs();
  }
  return n_geo_divs;
}

unsigned int CartogramInfo::n_insets() const
{
  return inset_states_.size();
}

std::map<std::string, InsetState> &CartogramInfo::ref_to_inset_states()
{
  return inset_states_;
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

  const double mean_density =
    total_target_area_with_data / total_start_area_with_data;

  std::cerr << "Mean density: " << mean_density << std::endl;

  // Calculate threshold for small target areas. For GeoDivs below this
  // threshold, the target area is scaled up for easier calculation.
  const double small_target_area_threshold =
    total_target_area_with_data * small_area_threshold_frac;

  std::cerr << "Using Small target area threshold: "
            << small_target_area_threshold << std::endl;

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
    // small_target_area_threshold if not all target areas are initially
    // missing or zero
    if (small_target_area_threshold > 0.0) {
      std::cerr << "Replacing small target areas..." << std::endl;
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
        if (
          (target_area >= 0.0) &&
          (target_area <= small_target_area_threshold)) {

          // Do not allow the replacement target area to be smaller than the
          // GeoDiv's target area
          double gd_specific_replacement_target_area = std::max(
            std::min(replacement_target_area, gd.area() * mean_density),
            target_area);
          inset_state.replace_target_area(
            gd.id(),
            gd_specific_replacement_target_area);
          std::cerr << gd.id() << ": " << target_area << " to "
                    << gd_specific_replacement_target_area
                    << " Area: " << gd.area() << "\n";

          // Update total target area
          total_target_area_with_data +=
            (gd_specific_replacement_target_area - target_area);
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
            const double adjusted_mean_density =
              total_target_area_with_data / total_start_area_with_data;
            new_target_area = adjusted_mean_density * gd.area();
          }
          inset_state.replace_target_area(gd.id(), new_target_area);
        }
      }
    }
  }
}

std::string CartogramInfo::set_map_name(const std::string &map_name)
{
  map_name_ = map_name;
  if (map_name_.find_last_of("/\\") != std::string::npos) {
    map_name_ = map_name_.substr(map_name_.find_last_of("/\\") + 1);
  }
  if (map_name_.find('.') != std::string::npos) {
    map_name_ = map_name_.substr(0, map_name_.find('.'));
  }
  return map_name_;
}

void CartogramInfo::write_csv(const std::string &csv_file_name) {
  // Write a csv file with the current target areas
  std::ofstream out_file_csv;
  out_file_csv.open(csv_file_name + ".csv");
  if (!out_file_csv) {
    std::cerr
      << "ERROR writing CSV: failed to open " << csv_file_name << ".csv"
      << std::endl;
  }

  // Each vector of strings will represent one row, starting with column names
  std::vector<std::vector<std::string> > csv_rows(1);

  csv_rows[0].push_back(id_header_);
  csv_rows[0].push_back("Target Area");

  // Fill up the rows with the IDs and target areas
  for (const auto &[id, inset_pos] : gd_to_inset_) {
    const auto &inset_state = inset_states_.at(inset_pos);
    const auto target_area = inset_state.target_area_at(id);
    csv_rows.push_back({id, std::to_string(target_area)});
  }

  // Write to CSV object
  auto writer = csv::make_csv_writer(out_file_csv);
  for (const auto &row : csv_rows) {
    writer << row;
  }

  // Close out_file and exit
  out_file_csv.close();
}