#ifndef CARTOGRAM_INFO_HPP_
#define CARTOGRAM_INFO_HPP_

#include "argparse.hpp"
#include "inset_state.hpp"

class CartogramInfo
{
private:
  std::map<std::string, std::string> gd_to_inset_;
  std::string id_header_;
  std::set<std::string> ids_in_visual_variables_file_;
  std::map<std::string, InsetState> inset_states_;
  bool is_world_map_;
  std::string map_name_;

  // TODO: We assume that either all external rings are counterclockwise or
  //       all are clockwise. This dichotomy covers most geospatial boundary
  //       files in the wild, but it would still be sensible to allow cases
  //       where there are external rings with opposite winding directions.
  bool original_ext_ring_is_clockwise_{};
  std::string visual_variable_file_;
  nlohmann::json cgal_to_json(bool = false);

public:
  explicit CartogramInfo(bool, const std::string &);
  [[nodiscard]] double cart_initial_total_target_area() const;
  [[nodiscard]] double area() const;
  [[nodiscard]] bool is_world_map() const;
  void json_to_geojson(
    const nlohmann::json &,
    nlohmann::ordered_json &,
    const nlohmann::json &);
  [[nodiscard]] unsigned int n_geo_divs() const;
  [[nodiscard]] unsigned int n_insets() const;
  void read_csv(const argparse::ArgumentParser &);
  void read_geojson(const std::string &, bool, std::string &);
  std::map<std::string, InsetState> &ref_to_inset_states();
  void replace_missing_and_zero_target_areas();
  std::string set_map_name(const std::string &);
  void shift_insets_to_target_position(bool = false);
  void write_csv(const std::string &csv_file_name);
  void write_geojson(const std::string &, const std::string &, bool = false);
};

#endif // CARTOGRAM_INFO_HPP_
