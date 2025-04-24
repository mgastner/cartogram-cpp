#include "inset_state.hpp"
#include "write_image.hpp"

// ======================== Grid Cell ========================

// Returns edge points of a grid cell and treats
// them as a polygon.
Polygon InsetState::grid_cell_edge_points(
  unsigned int x,
  unsigned int y,
  unsigned int cell_width,
  bool plot_equal_area_map) const
{
  Polygon cell_edge_points;
  const boost::multi_array<Point, 2> &proj =
    plot_equal_area_map ? identity_proj_ : cum_proj_;

  // Horizontal lower edge points
  for (unsigned int i = x; i < x + cell_width; ++i) {
    cell_edge_points.push_back(proj[i][y]);
  }

  // Vertical right edge points
  for (unsigned int i = y; i < y + cell_width; ++i) {
    cell_edge_points.push_back(proj[x + cell_width][i]);
  }

  // Horizontal upper edge points
  for (unsigned int i = x + cell_width; i > x; --i) {
    cell_edge_points.push_back(proj[i][y + cell_width]);
  }

  // Vertical left edge points
  for (unsigned int i = y + cell_width; i > y; --i) {
    cell_edge_points.push_back(proj[x][i]);
  }

  // Complete the polygon by making the first and last point the same
  cell_edge_points.push_back(cell_edge_points[0]);

  return cell_edge_points;
}

// Returns grid cell area based on edge points
double grid_cell_area(
  unsigned int x,
  unsigned int y,
  unsigned int cell_width,
  InsetState &inset_state)
{
  // Taking absolule to ensure we get the area irrespective of direction
  return abs(inset_state.grid_cell_edge_points(x, y, cell_width).area());
}

// Returns the largest and smallest grid cell area to be used for
// grid heatmap generation
std::pair<double, double> max_and_min_grid_cell_area(
  unsigned int cell_width,
  InsetState &inset_state)
{
  // Initialize max and min area
  double max_area = -dbl_inf;
  double min_area = dbl_inf;

  // Iterate over grid cells
  for (unsigned int i = 0; i < inset_state.lx() - cell_width;
       i += cell_width) {
    for (unsigned int j = 0; j < inset_state.ly() - cell_width;
         j += cell_width) {
      const double area = grid_cell_area(i, j, cell_width, inset_state);
      max_area = std::max(max_area, area);
      min_area = std::min(min_area, area);
    }
  }
  return std::make_pair(max_area, min_area);
}

// Returns the index of largest and smallest grid cell area to be used for
// grid heatmap generation
std::pair<Point, Point> max_and_min_grid_cell_area_index(
  unsigned int cell_width,
  InsetState &inset_state)
{
  // Initialize max and min area
  double max_area = -dbl_inf;
  double min_area = dbl_inf;
  unsigned int max_i = 0;
  unsigned int max_j = 0;
  unsigned int min_i = 0;
  unsigned int min_j = 0;

  // Iterate over grid cells
  for (unsigned int i = 0; i < inset_state.lx() - cell_width;
       i += cell_width) {
    for (unsigned int j = 0; j < inset_state.ly() - cell_width;
         j += cell_width) {
      const double area = grid_cell_area(i, j, cell_width, inset_state);
      if (area > max_area) {
        max_area = area;
        max_i = i;
        max_j = j;
      }
      if (area < min_area) {
        min_area = area;
        min_i = i;
        min_j = j;
      }
    }
  }
  return std::make_pair(Point(max_i, max_j), Point(min_i, min_j));
}

double grid_cell_target_area(
  const unsigned int i,
  const unsigned int j,
  const double total_target_area,
  const double total_inset_area,
  InsetState &inset_state)
{
  const Polygon cell_edge_points =
    inset_state.grid_cell_edge_points(i, j, 1, false);
  const double cell_area = cell_edge_points.area();

  const double cell_target_area =
    (cell_area * total_target_area) / total_inset_area;

  return cell_target_area;
}

double grid_cell_target_area_per_km(
  const unsigned int i,
  const unsigned int j,
  const double total_target_area,
  const double total_inset_area,
  InsetState &inset_state)
{
  const double cell_target_area = grid_cell_target_area(
    i,
    j,
    total_target_area,
    total_inset_area,
    inset_state);
  const double cell_area_km = grid_cell_area_km(inset_state, i, j);

  const double cell_target_area_per_km = cell_target_area / cell_area_km;

  return cell_target_area_per_km;
}

double grid_cell_area_km(
  const InsetState &inset_state,
  const unsigned int i,
  const unsigned int j)
{
  const Polygon cell_edge_points =
    inset_state.grid_cell_edge_points(i, j, 1, true);
  const Polygon cell_edge_points_equal_area_projection =
    transform_to_equal_area_projection_coor(cell_edge_points, inset_state);
  const double cell_area = cell_edge_points_equal_area_projection.area();
  const double cell_area_km =
    equal_area_projection_area_to_earth_area(cell_area);

  return cell_area_km;
}

// Find and return the area per grid cell
double compute_per_grid_cell(const InsetState &inset_state)
{
  double per_grid_cell = 1;

  const std::vector<double> multipliers = {2, 5, 2};
  const double total_area = inset_state.total_inset_area();
  size_t multiplier_idx = 0;

  while (per_grid_cell < 0.015 * total_area) {
    per_grid_cell *= multipliers[multiplier_idx];
    multiplier_idx = (multiplier_idx + 1) % multipliers.size();
  }
  return per_grid_cell;
}