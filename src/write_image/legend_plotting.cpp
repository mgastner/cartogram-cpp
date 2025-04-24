#include "inset_state.hpp"
#include "write_image.hpp"

// ======================== Legend Plotting ========================

// Adds a square legend outside the map
void write_legend_on_surface(
  cairo_t *cr,
  bool equal_area_map,
  const InsetState &inset_state)
{
  const std::pair<unsigned int, unsigned int> legend_info =
    equal_area_map ? get_km_legend_length(inset_state)
                   : get_visual_variable_legend_length(inset_state);

  Point legend_pos(0, 0);
  const unsigned int grid_cell_value = legend_info.first;
  const unsigned int total_value = legend_info.second;
  const double grid_cell_length = sqrt(compute_per_grid_cell(inset_state));
  cairo_set_line_width(cr, 1);
  const Color color("black");
  cairo_set_source_rgb(cr, color.r, color.g, color.b);
  cairo_rectangle(
    cr,
    legend_pos.x(),
    legend_pos.y(),
    grid_cell_length,
    grid_cell_length);
  cairo_stroke(cr);

  // Determine position of text on map
  double x = legend_pos.x() + (grid_cell_length * 1.25);
  double grid_cell_y = legend_pos.y() + (grid_cell_length * 0.5);
  double total_y = grid_cell_y + (grid_cell_length * 0.75);

  // Adjust font size according to size of SVG
  unsigned int font_size = 8 * (inset_state.lx() / 256);
  cairo_select_font_face(
    cr,
    "Sans",
    CAIRO_FONT_SLANT_NORMAL,
    CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size(cr, font_size);

  const std::pair<std::string, std::string> legend_labels =
    get_legend_labels(grid_cell_value, total_value);
  std::string grid_cell_label = legend_labels.first;
  std::string total_label = legend_labels.second;

  // Add known units
  if (equal_area_map) {
    grid_cell_label += " km²";
    total_label += " km²";
  } else {
    grid_cell_label += " people";
    total_label += " people";
  }

  // Write grid value and total value text onto map
  cairo_move_to(cr, x, grid_cell_y + (font_size / 2.0));
  cairo_show_text(cr, grid_cell_label.c_str());

  cairo_move_to(cr, x, total_y + (font_size / 2.0));
  cairo_show_text(cr, total_label.c_str());
}

std::pair<unsigned int, unsigned int> get_km_legend_length(
  const InsetState &inset_state)
{
  const double cell_area_km = grid_cell_area_km(inset_state);
  const unsigned int grid_cell_area = get_nearest_nice_number_for_legend(
    cell_area_km * compute_per_grid_cell(inset_state));
  const unsigned int total_area =
    cell_area_km * inset_state.total_inset_area();

  return std::pair<unsigned int, unsigned int>(grid_cell_area, total_area);
}

std::pair<unsigned int, unsigned int> get_visual_variable_legend_length(
  const InsetState &inset_state)
{
  const unsigned int per_area =
    inset_state.initial_target_area() / inset_state.total_inset_area();
  const unsigned int grid_cell_area = get_nearest_nice_number_for_legend(
    per_area * compute_per_grid_cell(inset_state));
  const unsigned int total_area = inset_state.initial_target_area();

  return std::pair<unsigned int, unsigned int>(grid_cell_area, total_area);
}

// Generate text labels for use in legend grid display
std::pair<std::string, std::string> get_legend_labels(
  unsigned int grid_cell_value,
  unsigned int total_value)
{
  // Display value per grid cell in billions/millions/thousands
  std::string grid_cell_label = "";
  if (grid_cell_value >= 1000000000) {
    const int billions = grid_cell_value / 1000000000;
    grid_cell_label = std::to_string(billions) + "B";
  } else if (grid_cell_value >= 1000000) {
    const int millions = grid_cell_value / 1000000;
    grid_cell_label = std::to_string(millions) + "M";
  } else if (grid_cell_value >= 1000) {
    const int thousands = grid_cell_value / 1000;
    grid_cell_label = std::to_string(thousands) + "K";
  } else {
    grid_cell_label = std::to_string(grid_cell_value);
  }

  // Display total value in billions/millions/thousands with 1 decimal place
  std::string total_label = "Total: ";
  std::stringstream sstream;
  if (total_value >= 1000000000) {
    const double billions = (double)total_value / 1000000000;
    sstream << std::fixed << std::setprecision(1) << billions;
    total_label += sstream.str() + "B";
  } else if (total_value >= 1000000) {
    const double millions = (double)total_value / 1000000;
    sstream << std::fixed << std::setprecision(1) << millions;
    total_label += sstream.str() + "M";
  } else if (total_value >= 1000) {
    const double thousands = (double)total_value / 1000;
    sstream << std::fixed << std::setprecision(1) << thousands;
    total_label += sstream.str() + "K";
  } else {
    total_label += std::to_string(total_value);
  }

  return std::pair<std::string, std::string>(grid_cell_label, total_label);
}

// Find the nearest matching nice number and value for use in legend
int get_nearest_nice_number_for_legend(int value)
{
  const int NICE_NUMBERS[4] = {1, 2, 5, 10};

  int new_value = 99;
  const int scale = std::floor(std::log10(value));
  const double value_first_digit =
    value / pow(10.0, scale);  // Get first digit of value in decimals
  double value_diff = abs(value_first_digit - new_value);

  // Loop through array of nice numbers to find the closest matching one
  for (int num : NICE_NUMBERS) {
    if (abs(value_first_digit - num) < value_diff) {
      value_diff = abs(value_first_digit - num);
      new_value = num;
    }
  }

  // Get the real nice number by multiplying the gotten value with the scale
  new_value = new_value * pow(10.0, scale);
  return new_value;
}

// CAN REMOVE
// Given coordinates in lx by ly coordinate system, returns the corresponding
// coordinates in the equal_area_projection projection coordinate system
Polygon transform_to_equal_area_projection_coor(
  Polygon edge_points,
  const InsetState &inset_state)
{
  Transformation scale(CGAL::SCALING, inset_state.latt_const());

  Polygon cell_edge_points_equal_area_projection =
    transform(scale, edge_points);
  return cell_edge_points_equal_area_projection;
}

// Given area in the equal_area_projection projection coordinate system,
// returns the corresponding area in the square km^2
double equal_area_projection_area_to_earth_area(
  const double equal_area_projection_area)
{
  return (equal_area_projection_area * earth_surface_area) / (4 * pi);
}

// Given area in square km^2 on Earth's surface,
// returns the corresponding area in the equal_area_projection projection
// coordinate system
double earth_area_to_equal_area_projection_area(const double earth_area)
{
  return (earth_area * 4 * pi) / earth_surface_area;
}