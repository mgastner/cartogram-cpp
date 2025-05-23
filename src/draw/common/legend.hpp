#pragma once

#include "basic_figures.hpp"
#include "canvas.hpp"
#include "color.hpp"
#include "geometry.hpp"
#include "inset_state.hpp"

// ======================== Legend Plotting ========================

// Find the nearest matching nice number and value for use in legend
static unsigned long get_nearest_nice_number_for_legend(double value)
{
  const unsigned long NICE_NUMBERS[4] = {1, 2, 5, 10};

  unsigned long new_value = 99;
  const int scale = static_cast<int>(std::floor(std::log10(value)));
  const double value_first_digit =
    value / pow(10.0, scale);  // Get first digit of value in decimals
  double value_diff = abs(value_first_digit - static_cast<double>(new_value));

  // Loop through array of nice numbers to find the closest matching one
  for (auto num : NICE_NUMBERS) {
    if (abs(value_first_digit - static_cast<double>(num)) < value_diff) {
      value_diff = abs(value_first_digit - static_cast<double>(num));
      new_value = num;
    }
  }

  // Get the real nice number by multiplying the gotten value with the scale
  new_value *= static_cast<unsigned long>(std::pow(10.0, scale) + 0.5);
  return new_value;
}

static std::pair<unsigned long, unsigned long> get_km_legend_length(
  const InsetState &inset_state)
{
  const double cell_area_km = grid_cell_area_km(inset_state);
  const unsigned long grid_cell_area = get_nearest_nice_number_for_legend(
    cell_area_km * compute_per_grid_cell(inset_state));
  const unsigned long total_area =
    static_cast<unsigned long>(cell_area_km * inset_state.total_inset_area());

  return std::pair<unsigned long, unsigned long>(grid_cell_area, total_area);
}

static std::pair<unsigned long, unsigned long> get_visual_variable_legend_length(
  const InsetState &inset_state)
{
  const double per_area =
    inset_state.initial_target_area() / inset_state.total_inset_area();
  const unsigned long grid_cell_area = get_nearest_nice_number_for_legend(
    per_area * compute_per_grid_cell(inset_state));
  const unsigned long total_area = static_cast<unsigned long>(inset_state.initial_target_area());

  return std::pair<unsigned long, unsigned long>(grid_cell_area, total_area);
}

// Generate text labels for use in legend grid display
static std::pair<std::string, std::string> get_legend_labels(
  unsigned long grid_cell_value,
  unsigned long total_value)
{
  // Display value per grid cell in billions/millions/thousands
  std::string grid_cell_label = "";
  if (grid_cell_value >= 1000000000UL) {
    const unsigned long billions = grid_cell_value / 1000000000UL;
    grid_cell_label = std::to_string(billions) + "B";
  } else if (grid_cell_value >= 1000000UL) {
    const unsigned long millions = grid_cell_value / 1000000UL;
    grid_cell_label = std::to_string(millions) + "M";
  } else if (grid_cell_value >= 1000UL) {
    const unsigned long thousands = grid_cell_value / 1000UL;
    grid_cell_label = std::to_string(thousands) + "K";
  } else {
    grid_cell_label = std::to_string(grid_cell_value);
  }

  // Display total value in billions/millions/thousands with 1 decimal place
  std::string total_label = "Total: ";
  std::stringstream sstream;
  if (total_value >= 1000000000) {
    const double billions = static_cast<double>(total_value) / 1000000000;
    sstream << std::fixed << std::setprecision(1) << billions;
    total_label += sstream.str() + "B";
  } else if (total_value >= 1000000) {
    const double millions = static_cast<double>(total_value) / 1000000;
    sstream << std::fixed << std::setprecision(1) << millions;
    total_label += sstream.str() + "M";
  } else if (total_value >= 1000) {
    const double thousands = static_cast<double>(total_value) / 1000;
    sstream << std::fixed << std::setprecision(1) << thousands;
    total_label += sstream.str() + "K";
  } else {
    total_label += std::to_string(total_value);
  }

  return std::pair<std::string, std::string>(grid_cell_label, total_label);
}

static void write_legend(
  Canvas &cvs,
  bool equal_area_map,
  const InsetState &inset_state)
{
  auto [grid_cell_val, total_val] =
    equal_area_map ? get_km_legend_length(inset_state)
                   : get_visual_variable_legend_length(inset_state);

  const Point legend_pos{0, 0};
  const double cell_len = std::sqrt(compute_per_grid_cell(inset_state));

  cvs.set_stroke(Color{"#000000"}, 1.0);
  cvs.rectangle(
    legend_pos.x(),
    legend_pos.y(),
    cell_len,
    cell_len);

  double x_text = legend_pos.x() + cell_len * 1.25;
  double y_grid_label = legend_pos.y() + cell_len * 0.50;
  double y_total_label = y_grid_label + cell_len * 0.75;

  unsigned font_px = 8 * (inset_state.lx() / 256);

  auto labels = get_legend_labels(grid_cell_val, total_val);
  std::string grid_lbl = labels.first;
  std::string total_lbl = labels.second;

  grid_lbl += equal_area_map ? " km²" : " people";
  total_lbl += equal_area_map ? " km²" : " people";

  cvs.text(x_text, y_grid_label + font_px * 0.5, grid_lbl, font_px);
  cvs.text(x_text, y_total_label + font_px * 0.5, total_lbl, font_px);
}
