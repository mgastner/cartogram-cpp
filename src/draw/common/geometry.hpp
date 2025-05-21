#pragma once

#include "basic_figures.hpp"
#include "inset_state.hpp"

// TODO: IS THERE A CGAL WAY OF DETERMINING WHETHER THE LABEL'S BOUNDING
//       BOX IS COMPLETELY CONTAINED IN THE POLYGON?

// TODO: SHOULD THE CRITERION FOR PRINTING A LABEL BE THAT IT FITS INSIDE THE
//       POLYGON WITH HOLES? THAT CRITERION WOULD BE MORE RESTRICTIVE THAN
//       FITTING INSIDE THE EXTERIOR RING.
static bool all_points_inside_exterior_ring(
  const std::vector<Point> &pts,
  const Polygon_with_holes &pwh)
{
  for (const auto &pt : pts) {
    if (pwh.outer_boundary().has_on_unbounded_side(pt)) {
      return false;
    }
  }
  return true;
}

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

// CAN REMOVE
// Given coordinates in lx by ly coordinate system, returns the corresponding
// coordinates in the equal_area_projection projection coordinate system
static Polygon transform_to_equal_area_projection_coor(
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
static double equal_area_projection_area_to_earth_area(
  const double equal_area_projection_area)
{
  return (equal_area_projection_area * earth_surface_area) / (4 * pi);
}

static double grid_cell_area_km(
  const InsetState &inset_state,
  const unsigned int i = 0,
  const unsigned int j = 0)
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
static double compute_per_grid_cell(const InsetState &inset_state)
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
