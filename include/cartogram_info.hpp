#ifndef CARTOGRAM_INFO_HPP_
#define CARTOGRAM_INFO_HPP_

#include "inset_state.hpp"
#include "parse_arguments.hpp"
#include "time_tracker.hpp"

class CartogramInfo
{
private:
  std::map<std::string, std::string> gd_to_inset_;
  Arguments args_;
  std::string id_header_;
  std::set<std::string> ids_in_visual_variables_file_;
  std::vector<std::string> initial_id_order_;
  std::vector<InsetState> inset_states_;
  bool is_world_map_;
  std::string map_name_;
  std::map<std::string, std::vector<std::string>> unique_properties_map_;
  int id_col_;

  // TODO: We assume that either all external rings are counterclockwise or
  //       all are clockwise. This dichotomy covers most geospatial boundary
  //       files in the wild, but it would still be sensible to allow cases
  //       where there are external rings with opposite winding directions.
  bool original_ext_ring_is_clockwise_{};
  nlohmann::json cgal_to_json(bool = false);

  std::string crs_;

  TimeTracker timer;

  // Make default constructor private so that only
  // CartogramInfo(const std::string, Arguments) can be called as constructor
  CartogramInfo();

public:
  explicit CartogramInfo(Arguments);
  [[nodiscard]] double cart_initial_total_target_area() const;
  void construct_inset_state_from_geodivs(const simdjson::dom::element &);
  bool converged() const;
  [[nodiscard]] double area() const;
  [[nodiscard]] bool is_world_map() const;
  void json_to_geojson(
    const nlohmann::json &,
    nlohmann::ordered_json &,
    const nlohmann::json &);
  [[nodiscard]] size_t n_geo_divs() const;
  [[nodiscard]] size_t n_insets() const;

  void read_csv();
  void read_geojson();
  std::vector<InsetState> &ref_to_inset_states();
  void relocate_geodivs_based_on_inset_pos(
    const std::map<std::string, std::map<std::string, std::string>> &);
  void replace_missing_and_zero_target_areas();

  // Rescale insets in correct proportion to each other
  void rescale_insets();

  std::string set_map_name(const std::string &);
  void reposition_insets(bool output_to_stdout = false);

  void plot_input();
  void preprocess();
  void print_time_report();
  void project_to_equal_area();

  std::string match_id_columns(const std::optional<std::string> &);
  void update_id_header_info(const std::string &);
  void write_csv(const std::string &csv_file_name);
  void write_geojson(const std::string &);
  void write_shifted_insets();
  InsetState convert_to_inset_state();
  void write_svg(const std::string &suffix = "");
};

#endif  // CARTOGRAM_INFO_HPP_
