#include "inset_state.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <queue>
#include <set>
#include <vector>

using Coordinate = std::pair<int, int>;
using GridCoordinates = std::vector<Coordinate>;

// ---------------------------------------------------------------------
// Polygon Clipping (Sutherland–Hodgman)
// ---------------------------------------------------------------------

// Clips a polygon using a generic Sutherland–Hodgman routine
static Polygon clip_polygon_sutherland_hodgman(
  const Polygon &polygon,
  const std::function<bool(const Point &)> &is_inside,
  const std::function<Point(const Point &, const Point &)>
    &compute_intersection)
{
  if (polygon.size() == 0)
    return polygon;

  Polygon clipped;

  // Start with the last vertex to close the polygon
  Point prev_point = polygon[polygon.size() - 1];
  bool prev_inside = is_inside(prev_point);

  for (unsigned int i = 0; i < polygon.size(); ++i) {
    Point curr_point = polygon[i];
    bool curr_inside = is_inside(curr_point);

    if (curr_inside) {
      if (!prev_inside)
        clipped.push_back(compute_intersection(prev_point, curr_point));
      clipped.push_back(curr_point);
    } else if (prev_inside) {
      clipped.push_back(compute_intersection(prev_point, curr_point));
    }

    prev_point = curr_point;
    prev_inside = curr_inside;
  }
  return clipped;
}

// Clip polygon by a vertical line x = value
// If left_side is true, keeps points with x >= value; otherwise, x <= value
Polygon clip_polygon_by_vertical_line(
  const Polygon &polygon,
  double x,
  bool left_side)
{
  auto is_inside = [x, left_side](const Point &p) -> bool {
    return left_side ? (p.x() >= x) : (p.x() <= x);
  };
  auto compute_intersection = [x](const Point &p, const Point &q) -> Point {
    double y = p.y() + (q.y() - p.y()) * (x - p.x()) / (q.x() - p.x());
    return Point(x, y);
  };
  return clip_polygon_sutherland_hodgman(
    polygon,
    is_inside,
    compute_intersection);
}

// Clip polygon by a horizontal line y = value
// If bottom_side is true, keeps points with y >= value; otherwise, y <= value
Polygon clip_polygon_by_horizontal_line(
  const Polygon &polygon,
  double y,
  bool bottom_side)
{
  auto is_inside = [y, bottom_side](const Point &p) -> bool {
    return bottom_side ? (p.y() >= y) : (p.y() <= y);
  };
  auto compute_intersection = [y](const Point &p, const Point &q) -> Point {
    double x = p.x() + (q.x() - p.x()) * (y - p.y()) / (q.y() - p.y());
    return Point(x, y);
  };
  return clip_polygon_sutherland_hodgman(
    polygon,
    is_inside,
    compute_intersection);
}

// Compute the overlapping area between a square and a polygon
double compute_square_polygon_overlap_area(
  const Polygon &square,
  const Polygon &polygon)
{
  const Bbox bbox = square.bbox();
  Polygon clipped = polygon;
  clipped = clip_polygon_by_vertical_line(clipped, bbox.xmin(), true);
  clipped = clip_polygon_by_vertical_line(clipped, bbox.xmax(), false);
  clipped = clip_polygon_by_horizontal_line(clipped, bbox.ymin(), true);
  clipped = clip_polygon_by_horizontal_line(clipped, bbox.ymax(), false);
  return std::abs(clipped.area());
}

// Compute the overlapping area between a square and a polygon with holes
double compute_square_polygon_with_holes_overlap_area(
  const Polygon &square,
  const Polygon_with_holes &pwh)
{
  double area = 0.0;
  area += compute_square_polygon_overlap_area(square, pwh.outer_boundary());
  for (auto hole_it = pwh.holes_begin(); hole_it != pwh.holes_end();
       ++hole_it) {
    area -= compute_square_polygon_overlap_area(square, *hole_it);
  }
  return area;
}

// ---------------------------------------------------------------------
// Supercover Line and Edge Rasterization
// ---------------------------------------------------------------------

// Compute a supercover line (returns grid cells traversed by the
// line) based on Amanatides–Woo's line algorithm. Perfectly accurate
GridCoordinates compute_supercover_line(
  double x0,
  double y0,
  double x1,
  double y1)
{
  GridCoordinates cells;

  // Determine the starting and ending grid cell coordinates
  int cell_x = static_cast<int>(std::floor(x0));
  int cell_y = static_cast<int>(std::floor(y0));
  int end_cell_x = static_cast<int>(std::floor(x1));
  int end_cell_y = static_cast<int>(std::floor(y1));

  cells.push_back({cell_x, cell_y});

  // Ray direction
  double delta_x = x1 - x0;
  double delta_y = y1 - y0;

  // Step direction
  int step_x = (delta_x > 0) ? 1 : (delta_x < 0 ? -1 : 0);
  int step_y = (delta_y > 0) ? 1 : (delta_y < 0 ? -1 : 0);

  // These represent how far along the ray we must move (in t, where 0 <= t <=
  // 1) to cross a cell boundary in the respective direction
  double t_delta_x = (delta_x != 0) ? 1.0 / std::abs(delta_x)
                                    : std::numeric_limits<double>::infinity();
  double t_delta_y = (delta_y != 0) ? 1.0 / std::abs(delta_y)
                                    : std::numeric_limits<double>::infinity();

  // The parametric distance along the ray until we reach the first vertical
  // (or horizontal) grid line
  double t_max_x, t_max_y;
  if (delta_x != 0) {
    if (step_x > 0)
      t_max_x = ((cell_x + 1) - x0) / delta_x;
    else
      t_max_x = (x0 - cell_x) / -delta_x;
  } else {
    t_max_x = std::numeric_limits<double>::infinity();
  }

  if (delta_y != 0) {
    if (step_y > 0)
      t_max_y = ((cell_y + 1) - y0) / delta_y;
    else
      t_max_y = (y0 - cell_y) / -delta_y;
  } else {
    t_max_y = std::numeric_limits<double>::infinity();
  }

  double t = 0.0;
  // Traverse the grid until t exceeds 1 (the end of the segment) or we hit the
  // end cell
  while (t <= 1.0) {
    if (t_max_x < t_max_y) {
      cell_x += step_x;
      t = t_max_x;
      t_max_x += t_delta_x;
    } else {
      cell_y += step_y;
      t = t_max_y;
      t_max_y += t_delta_y;
    }
    if (t > 1.0)
      break;
    cells.push_back({cell_x, cell_y});
    if (cell_x == end_cell_x && cell_y == end_cell_y)
      break;
  }

  return cells;
}

// Runs the supercover line algorithm on the edges of a polygon
GridCoordinates rasterize_polygon_edges(const Polygon &polygon)
{
  GridCoordinates cells;
  if (polygon.size() < 2)
    return cells;
  for (std::size_t i = 0; i < polygon.size(); ++i) {
    Point p0 = polygon[i];
    Point p1 = polygon[(i + 1) % polygon.size()];
    GridCoordinates line_cells =
      compute_supercover_line(p0.x(), p0.y(), p1.x(), p1.y());
    cells.insert(cells.end(), line_cells.begin(), line_cells.end());
  }
  return cells;
}

// ---------------------------------------------------------------------
// Point-in-Polygon Check for Polygon_with_holes
// ---------------------------------------------------------------------

bool point_inside_polygon(const Polygon &polygon, const Point &pt)
{
  return polygon.bounded_side(pt) == CGAL::ON_BOUNDED_SIDE;
}

bool point_inside_polygon_with_holes(
  const Polygon_with_holes &pwh,
  const Point &pt)
{
  if (!point_inside_polygon(pwh.outer_boundary(), pt)) {
    return false;
  }
  for (auto hole_it = pwh.holes_begin(); hole_it != pwh.holes_end();
       ++hole_it) {
    if (point_inside_polygon(*hole_it, pt)) {
      return false;
    }
  }
  return true;
}

// ---------------------------------------------------------------------
// Data Structures for Polygon Information
// ---------------------------------------------------------------------

struct PolygonInfo {
  unsigned int gd_id;
  unsigned int pwh_tot_id;
  unsigned int pwh_id;
  bool is_hole;
  unsigned int hole_id;
};

// Information about which Polygon lies in Cell (x, y)
struct GridCellPolygon {
  int x;
  int y;
  PolygonInfo poly_info;
};

// ---------------------------------------------------------------------
// Helper Functions for fill_with_density_clip
// ---------------------------------------------------------------------

// Runs the supercover line algorithm on the edges of the Map
// * For each cell, stores which pwhs' edges are present in the cell
// * Whether the cell is an edge cell
// * Stores metadata for each pwh
void process_geo_divisions_edge_info(
  boost::multi_array<bool, 2> &is_edge,
  boost::multi_array<std::vector<PolygonInfo>, 2> &grid_cell_polygons,
  std::vector<PolygonInfo> &all_pwh_info,
  InsetState &inset_state)
{
  std::vector<GridCellPolygon> global_edge_info;
  global_edge_info.reserve(10000);

#pragma omp parallel
  {
    std::vector<GridCellPolygon> local_edge_info;
#pragma omp for nowait schedule(dynamic)
    for (unsigned int gd_id = 0; gd_id < inset_state.geo_divs().size();
         ++gd_id) {
      const GeoDiv &gd = inset_state.geo_divs()[gd_id];
      for (unsigned int pwh_id = 0; pwh_id < gd.n_polygons_with_holes();
           ++pwh_id) {
        const Polygon_with_holes &pwh = gd.polygons_with_holes()[pwh_id];
        const PolygonInfo outer_info{
          gd_id,
          static_cast<unsigned int>(all_pwh_info.size()),
          pwh_id,
          false,
          0};
        const GridCoordinates outer_cells =
          rasterize_polygon_edges(pwh.outer_boundary());
        for (const auto &cell : outer_cells) {
          const auto &[x, y] = cell;
          local_edge_info.push_back({x, y, outer_info});
        }
        for (unsigned int hole_id = 0; hole_id < pwh.number_of_holes();
             ++hole_id) {
          const Polygon &hole = pwh.holes()[hole_id];
          const GridCoordinates hole_cells = rasterize_polygon_edges(hole);
          const PolygonInfo hole_info{
            gd_id,
            static_cast<unsigned int>(all_pwh_info.size()),
            pwh_id,
            true,
            hole_id};
          for (const auto &cell : hole_cells) {
            const auto &[x, y] = cell;
            local_edge_info.push_back({x, y, hole_info});
          }
        }
#pragma omp critical
        {
          all_pwh_info.push_back(outer_info);
        }
      }
    }
#pragma omp critical
    {
      global_edge_info.insert(
        global_edge_info.end(),
        local_edge_info.begin(),
        local_edge_info.end());
    }
  }

  for (const auto &gcp : global_edge_info) {
    const auto &[x, y, poly_info] = gcp;
    is_edge[x][y] = true;
    grid_cell_polygons[x][y].push_back(poly_info);
  }
}

// For all non-edge cells, compute connected components
// and assign a unique ID to each connected component
// This is useful later to classify all cells of a connected
// component to a single pwh or to ocean
void compute_connected_components(
  boost::multi_array<int, 2> &comp,
  int &current_comp_id,
  const boost::multi_array<bool, 2> &is_edge,
  InsetState &inset_state)
{
  current_comp_id = 0;
  int dx[4] = {0, 0, 1, -1};
  int dy[4] = {1, -1, 0, 0};
  for (int x = 0; x < static_cast<int>(inset_state.lx()); ++x) {
    for (int y = 0; y < static_cast<int>(inset_state.ly()); ++y) {
      if (!is_edge[x][y] && comp[x][y] == -1) {
        std::queue<Coordinate> q;
        q.push({x, y});
        comp[x][y] = current_comp_id;
        while (!q.empty()) {
          auto [cx, cy] = q.front();
          q.pop();
          for (int d = 0; d < 4; ++d) {
            int nx = cx + dx[d];
            int ny = cy + dy[d];
            if (
              nx < 0 || nx >= static_cast<int>(inset_state.lx()) || ny < 0 ||
              ny >= static_cast<int>(inset_state.ly()))
              continue;
            if (!is_edge[nx][ny] && comp[nx][ny] == -1) {
              comp[nx][ny] = current_comp_id;
              q.push({nx, ny});
            }
          }
        }
        ++current_comp_id;
      }
    }
  }
}

// If inside, for each connected component (cc), maps it to corresponding pwh
// For an unmapped cell of a cc, if the neighboring cell is a edge cell,
// then we find the pwh of that edge, and run point-in-polygon check
// to see if the cc is inside the pwh. If so, we map the cc to that pwh
void map_components_to_polygon_ids(
  std::vector<bool> &is_comp_used,
  std::vector<int> &comp_id_to_pwh_id,
  const boost::multi_array<bool, 2> &is_edge,
  const boost::multi_array<std::vector<PolygonInfo>, 2> &grid_cell_polygons,
  const boost::multi_array<int, 2> &comp,
  const std::vector<PolygonInfo> &all_pwh_info,
  InsetState &inset_state)
{
  // Store the connected components that are outside the pwh when encountered
  // for the first time. Avoids running the expensive point-in-polygon check
  // multiple times
  std::vector<std::set<int>> is_outside(all_pwh_info.size());

  int dx[4] = {0, 0, 1, -1};
  int dy[4] = {1, -1, 0, 0};

  for (int x = 0; x < static_cast<int>(inset_state.lx()); ++x) {
    for (int y = 0; y < static_cast<int>(inset_state.ly()); ++y) {
      const int comp_id = comp[x][y];
      if (is_edge[x][y] || is_comp_used[comp_id]) {
        continue;
      }
      for (int d = 0; d < 4; ++d) {
        const int nx = x + dx[d];
        const int ny = y + dy[d];
        if (
          nx < 0 || nx >= static_cast<int>(inset_state.lx()) || ny < 0 ||
          ny >= static_cast<int>(inset_state.ly())) {
          continue;
        }
        if (is_edge[nx][ny]) {
          for (const auto &poly_info : grid_cell_polygons[nx][ny]) {
            const unsigned int pwh_tot_id = poly_info.pwh_tot_id;
            if (is_outside[pwh_tot_id].contains(comp_id)) {
              continue;
            }
            const unsigned int gd_id = poly_info.gd_id;
            const unsigned int pwh_id = poly_info.pwh_id;
            const bool is_hole = poly_info.is_hole;
            const auto &pwh =
              inset_state.geo_divs()[gd_id].polygons_with_holes()[pwh_id];
            const Point center = Point(x + 0.5, y + 0.5);
            if (is_hole) {
              const unsigned int &hole_id = poly_info.hole_id;
              auto &hole = pwh.holes()[hole_id];
              if (point_inside_polygon(hole, center)) {
                is_outside[pwh_tot_id].insert(comp_id);
              } else {
                is_comp_used[comp_id] = true;
                comp_id_to_pwh_id[comp_id] = pwh_tot_id;
              }
            } else {
              auto &outer_boundary = pwh.outer_boundary();
              if (point_inside_polygon(outer_boundary, center)) {
                is_comp_used[comp_id] = true;
                comp_id_to_pwh_id[comp_id] = pwh_tot_id;
              } else {
                is_outside[pwh_tot_id].insert(comp_id);
              }
            }
          }
        }
      }
    }
  }
}

// Computes area, using area computes the density metrics for each cell
// If a cell is fully inside a pwh, then we know the area in O(1)
// If a cell is an edge cell, then we need to run expensive clipping
// algorithm to compute the area (this is the bottleneck)
void compute_density_grid(
  boost::multi_array<double, 2> &area_filled,
  boost::multi_array<double, 2> &numer,
  boost::multi_array<double, 2> &denom,
  const boost::multi_array<bool, 2> &is_edge,
  const boost::multi_array<int, 2> &comp,
  const std::vector<bool> &is_comp_used,
  const std::vector<PolygonInfo> &all_pwh_info,
  const boost::multi_array<std::vector<PolygonInfo>, 2> &grid_cell_polygons,
  const std::vector<int> &comp_id_to_pwh_to_id,
  const std::vector<double> &gd_target_density,
  InsetState &inset_state)
{
  std::fill_n(area_filled.data(), area_filled.num_elements(), 0.0);
  std::fill_n(numer.data(), numer.num_elements(), 0.0);
  std::fill_n(denom.data(), denom.num_elements(), 0.0);

#pragma omp parallel for collapse(2) schedule(dynamic)
  for (int x = 0; x < static_cast<int>(inset_state.lx()); ++x) {
    for (int y = 0; y < static_cast<int>(inset_state.ly()); ++y) {
      const int comp_id = comp[x][y];
      // Skip ocean cells
      if (not is_edge[x][y] and not is_comp_used[comp_id]) {
        continue;
      }

      // Build cell polygon (square cell).
      Polygon cell;
      cell.push_back(Point(x, y));
      cell.push_back(Point(x + 1, y));
      cell.push_back(Point(x + 1, y + 1));
      cell.push_back(Point(x, y + 1));

      if (!is_edge[x][y] and is_comp_used[comp_id]) {
        // Fully inside a Polygon with holes
        const unsigned int pwh_tot_id = comp_id_to_pwh_to_id[comp_id];
        const PolygonInfo &poly_info = all_pwh_info[pwh_tot_id];
        const unsigned int gd_id = poly_info.gd_id;
        const std::string &gd_name = inset_state.geo_divs()[gd_id].id();
        const double weight = 1.0 * inset_state.area_error_at(gd_name);
#pragma omp atomic
        numer[x][y] += weight * gd_target_density[gd_id];
#pragma omp atomic
        denom[x][y] += weight;
        area_filled[x][y] = 1.0;
        continue;
      }

      const std::vector<PolygonInfo> &cell_poly_info =
        grid_cell_polygons[x][y];

      // Avoid double counting the same pwh
      std::set<int> visited_pwh_ids;
      for (const auto &poly_info : cell_poly_info) {
        const unsigned int pwh_tot_id = poly_info.pwh_tot_id;
        if (visited_pwh_ids.contains(pwh_tot_id)) {
          continue;
        }
        visited_pwh_ids.insert(pwh_tot_id);  // Mark as visited
        const unsigned int gd_id = poly_info.gd_id;
        const unsigned int pwh_id = poly_info.pwh_id;
        const std::string &gd_name = inset_state.geo_divs()[gd_id].id();
        const Polygon_with_holes &pwh =
          inset_state.geo_divs()[gd_id].polygons_with_holes()[pwh_id];

        // Bottleneck: Clip cell with polygon with holes
        const double intersect_area =
          compute_square_polygon_with_holes_overlap_area(cell, pwh);
#pragma omp atomic
        area_filled[x][y] += intersect_area;
        double weight = intersect_area * inset_state.area_error_at(gd_name);
#pragma omp atomic
        numer[x][y] += weight * gd_target_density[gd_id];
#pragma omp atomic
        denom[x][y] += weight;
      }
    }
  }
}

// ---------------------------------------------------------------------
// Testing and Debugging
// ---------------------------------------------------------------------

// Test function to compare the computed area with the actual area (brute
// force)
void test_areas(
  const boost::multi_array<double, 2> &area_found,
  InsetState &inset_state)
{
  const unsigned int lx = inset_state.lx();
  const unsigned int ly = inset_state.ly();
  boost::multi_array<double, 2> area_actual(boost::extents[lx][ly]);
  std::fill_n(area_actual.data(), area_actual.num_elements(), 0.0);
  for (unsigned int gd_id = 0; gd_id < inset_state.geo_divs().size();
       ++gd_id) {
    const GeoDiv &gd = inset_state.geo_divs()[gd_id];
    for (const auto &pwh : gd.polygons_with_holes()) {
      auto bbox = pwh.outer_boundary().bbox();
      for (unsigned int i = std::max(0, static_cast<int>(bbox.xmin()));
           i < std::min(lx, static_cast<unsigned int>(bbox.xmax()) + 1);
           ++i) {
        for (unsigned int j = std::max(0, static_cast<int>(bbox.ymin()));
             j < std::min(ly, static_cast<unsigned int>(bbox.ymax()) + 1);
             ++j) {

          // Create a 1x1 cell
          Polygon cell;
          cell.push_back(Point(i, j));
          cell.push_back(Point(i + 1, j));
          cell.push_back(Point(i + 1, j + 1));
          cell.push_back(Point(i, j + 1));

          double intersect_area_pwh =
            compute_square_polygon_with_holes_overlap_area(cell, pwh);
          area_actual[i][j] += intersect_area_pwh;
        }
      }
    }
  }
  // do m, mean absolute error, etc.
  double absolute_error = 0.0;
  double square_error = 0.0;
  const unsigned int num_cells = lx * ly;
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      absolute_error += std::abs(area_actual[i][j] - area_found[i][j]);
      square_error += std::pow(area_actual[i][j] - area_found[i][j], 2);
    }
  }
  const double mean_absolute_error = absolute_error / num_cells;
  const double root_mean_square_error = std::sqrt(square_error / num_cells);
  std::cout << "Num Integration: " << inset_state.n_finished_integrations()
            << std::endl;
  std::cout << "Mean Absolute Error: " << mean_absolute_error << std::endl;
  std::cout << "Root Mean Square Error: " << root_mean_square_error
            << std::endl
            << std::endl;
}

void InsetState::create_contiguity_graph() {

  // Create a copy of `inset_state` to avoid rescaling the original
  InsetState is_copy = *this;

  is_copy.rescale_map();

  boost::multi_array<bool, 2> is_edge(boost::extents[is_copy.lx()][is_copy.ly()]);
  std::fill_n(is_edge.data(), is_edge.num_elements(), false);
  boost::multi_array<std::vector<PolygonInfo>, 2> grid_cell_polygons(
    boost::extents[is_copy.lx()][is_copy.ly()]);
  std::vector<PolygonInfo> all_pwh_info;

  process_geo_divisions_edge_info(
    is_edge,
    grid_cell_polygons,
    all_pwh_info,
    is_copy);

  for (unsigned int x = 0; x < is_copy.lx(); ++x) {
    for (unsigned int y = 0; y < is_copy.ly(); ++y) {
      if (is_edge[x][y]) {

        // Iterate through
        std::vector<PolygonInfo> polygons_at_xy = grid_cell_polygons[x][y];

        for (size_t i = 0; i < polygons_at_xy.size(); ++i) {
          GeoDiv &gd_1 = geo_divs_[polygons_at_xy[i].gd_id];
          for (size_t j = i + 1; j < polygons_at_xy.size(); ++j) {
            GeoDiv &gd_2 = geo_divs_[polygons_at_xy[j].gd_id];

            // Update adjacency for both GeoDivs
            // GeoDiv.adjacent_to() already checks whether the GeoDivs are different
            gd_1.adjacent_to(gd_2.id());
            gd_2.adjacent_to(gd_1.id());
          }
        }

      }
    }
  }
}

void InsetState::fill_with_density_clip()
{
  timer.start("Fill with Density (Clipping Method)");

  // Step 1: Detect edges and store edge information
  boost::multi_array<bool, 2> is_edge(boost::extents[lx_][ly_]);
  std::fill_n(is_edge.data(), is_edge.num_elements(), false);
  boost::multi_array<std::vector<PolygonInfo>, 2> grid_cell_polygons(
    boost::extents[lx_][ly_]);
  std::vector<PolygonInfo> all_pwh_info;

  process_geo_divisions_edge_info(
    is_edge,
    grid_cell_polygons,
    all_pwh_info,
    *this);

  // Step 2: Compute connected components
  boost::multi_array<int, 2> comp(boost::extents[lx_][ly_]);
  std::fill_n(comp.data(), comp.num_elements(), -1);
  int current_comp_id = 0;
  compute_connected_components(comp, current_comp_id, is_edge, *this);

  // Step 3: Map components to polygon with holes IDs
  std::vector<bool> is_comp_used(current_comp_id, false);
  std::vector<int> comp_id_to_pwh_id(current_comp_id, -1);
  map_components_to_polygon_ids(
    is_comp_used,
    comp_id_to_pwh_id,
    is_edge,
    grid_cell_polygons,
    comp,
    all_pwh_info,
    *this);

  // Step 4: Compute density grid
  boost::multi_array<double, 2> area_filled(boost::extents[lx_][ly_]);
  boost::multi_array<double, 2> numer(boost::extents[lx_][ly_]);
  boost::multi_array<double, 2> denom(boost::extents[lx_][ly_]);

  std::vector<double> gd_target_density;
  gd_target_density.reserve(geo_divs_.size());
  for (const auto &gd : geo_divs_) {
    gd_target_density.push_back(target_area_at(gd.id()) / gd.area());
  }

#pragma omp parallel for collapse(2)
  for (unsigned int i = 0; i < lx_; ++i)
    for (unsigned int j = 0; j < ly_; ++j)
      rho_init_(i, j) = 0.0;

  compute_density_grid(
    area_filled,
    numer,
    denom,
    is_edge,
    comp,
    is_comp_used,
    all_pwh_info,
    grid_cell_polygons,
    comp_id_to_pwh_id,
    gd_target_density,
    *this);

  // test_areas(area_filled, *this);

  // Step 5: Compute ocean density and adjust remaining cells
  double ocean_density =
    (lx_ * ly_ - total_target_area()) / (lx_ * ly_ - total_inset_area());
  double ocean_area_error = std::abs(
    (lx_ * ly_ - total_inset_area()) / (lx_ * ly_ - total_target_area()) -
    1.0);
#pragma omp parallel for collapse(2)
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      double weight = (1.0 - area_filled[i][j]) * ocean_area_error;
#pragma omp atomic
      numer[i][j] += weight * ocean_density;
#pragma omp atomic
      denom[i][j] += weight;
    }
  }

#pragma omp parallel for collapse(2)
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      if (denom[i][j] > 0)
        rho_init_(i, j) = numer[i][j] / denom[i][j];
      else
        rho_init_(i, j) = ocean_density;
    }
  }

  // Step 6: Run FFTW to compute rho_ft_
  auto rho_begin = rho_init_.as_1d_array();
  auto rho_end = rho_begin + lx_ * ly_;
  auto [min_iter, max_iter] = std::minmax_element(rho_begin, rho_end);
  dens_min_ = *min_iter;
  dens_mean_ = ocean_density;
  dens_max_ = *max_iter;

  execute_fftw_fwd_plan();
  timer.stop("Fill with Density (Clipping Method)");
}
