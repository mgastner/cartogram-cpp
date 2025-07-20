#include "inset_state.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <queue>
#include <set>
#include <vector>

using Coordinate = std::pair<int, int>;
using EdgeEndpointIndices = std::pair<unsigned int, unsigned int>;
using GridCoordinates = std::vector<Coordinate>;
using GridCoordinatesWithEdgeEndpoints =
  std::vector<std::pair<Coordinate, EdgeEndpointIndices>>;

// ---------------------------------------------------------------------
// Clipping Functions
// ---------------------------------------------------------------------

// Open-path clipping using a Sutherland–Hodgman–like algorithm
// Here, 'path' is a polyline (open, not closed)
static Polygon clip_path_sutherland_hodgman(
  const Polygon &path,
  const std::function<bool(const Point &)> &is_inside,
  const std::function<Point(const Point &, const Point &)>
    &compute_intersection)
{
  if (path.is_empty()) {
    return path;
  }

  Polygon clipped;

  // Determine if the first point is inside
  bool prev_inside = is_inside(path[0]);
  Point prev_point = path[0];
  if (prev_inside) {
    clipped.push_back(prev_point);
  }

  // Process each consecutive pair
  for (unsigned int i = 1; i < path.size(); ++i) {
    Point curr_point = path[i];
    bool curr_inside = is_inside(curr_point);

    // If the segment crosses the clipping boundary:
    if (prev_inside != curr_inside) {

      Point inter_pt = compute_intersection(prev_point, curr_point);

      // If we're leaving the clipping region, add the intersection
      if (prev_inside && !curr_inside) {
        clipped.push_back(inter_pt);
      }

      // If we're entering the clipping region, add the intersection then the
      // point
      else if (!prev_inside && curr_inside) {
        clipped.push_back(inter_pt);
        clipped.push_back(curr_point);
      }
    }

    // If both endpoints are inside, simply add the current point
    else if (curr_inside) {
      clipped.push_back(curr_point);
    }

    prev_point = curr_point;
    prev_inside = curr_inside;
  }
  return clipped;
}

// Clip a path by a vertical line x = value
// If left_side is true, keep points with x >= value; otherwise, x <= value.
static Polygon clip_path_by_vertical_line(
  const Polygon &path,
  const double x,
  const bool left_side)
{
  auto is_inside = [x, left_side](const Point &p) -> bool {
    return left_side ? (p.x() >= x) : (p.x() <= x);
  };
  auto compute_intersection = [x](const Point &p, const Point &q) -> Point {
    double y = p.y() + (q.y() - p.y()) * (x - p.x()) / (q.x() - p.x());
    return Point(x, y);
  };
  return clip_path_sutherland_hodgman(path, is_inside, compute_intersection);
}

// Clip a path by a horizontal line y = value.
// If bottom_side is true, keep points with y >= value; otherwise, y <= value.
static Polygon clip_path_by_horizontal_line(
  const Polygon &path,
  const double y,
  const bool bottom_side)
{
  auto is_inside = [y, bottom_side](const Point &p) -> bool {
    return bottom_side ? (p.y() >= y) : (p.y() <= y);
  };
  auto compute_intersection = [y](const Point &p, const Point &q) -> Point {
    double x = p.x() + (q.x() - p.x()) * (y - p.y()) / (q.y() - p.y());
    return Point(x, y);
  };
  return clip_path_sutherland_hodgman(path, is_inside, compute_intersection);
}

static Polygon clip_path_by_rectangle(
  const Polygon &path,
  const Bbox &rect_bbox)
{
  Polygon clipped = path;
  clipped = clip_path_by_vertical_line(clipped, rect_bbox.xmin(), true);
  clipped = clip_path_by_vertical_line(clipped, rect_bbox.xmax(), false);
  clipped = clip_path_by_horizontal_line(clipped, rect_bbox.ymin(), true);
  clipped = clip_path_by_horizontal_line(clipped, rect_bbox.ymax(), false);
  return clipped;
}

// Similar to the idea used in `clip_path_sutherland_hodgman`, but for polygons
static Polygon clip_polygon_sutherland_hodgman(
  const Polygon &polygon,
  const std::function<bool(const Point &)> &is_inside,
  const std::function<Point(const Point &, const Point &)>
    &compute_intersection)
{
  if (polygon.is_empty())
    return polygon;

  Polygon clipped;

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

static Polygon clip_polygon_by_vertical_line(
  const Polygon &polygon,
  const double x,
  const bool left_side)
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

static Polygon clip_polygon_by_horizontal_line(
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

static Polygon clip_polygon_by_rectangle(
  const Polygon &polygon,
  const Bbox &rect_bbox)
{
  Polygon clipped = polygon;
  clipped = clip_polygon_by_vertical_line(clipped, rect_bbox.xmin(), true);
  clipped = clip_polygon_by_vertical_line(clipped, rect_bbox.xmax(), false);
  clipped = clip_polygon_by_horizontal_line(clipped, rect_bbox.ymin(), true);
  clipped = clip_polygon_by_horizontal_line(clipped, rect_bbox.ymax(), false);
  return clipped;
}

static double compute_polygon_rectangle_overlap_area(
  const Polygon &polygon,
  const Bbox &rect_bbox)
{
  Polygon clipped = clip_polygon_by_rectangle(polygon, rect_bbox);
  return std::abs(clipped.area());
}

// Compute the overlapping area between a square and a polygon with holes
static double compute_pwh_rectangle_overlap_area(
  const Polygon_with_holes &pwh,
  const Bbox &rect_bbox)
{
  double area = 0.0;
  area +=
    compute_polygon_rectangle_overlap_area(pwh.outer_boundary(), rect_bbox);
  for (auto hole_it = pwh.holes_begin(); hole_it != pwh.holes_end();
       ++hole_it) {
    area -= compute_polygon_rectangle_overlap_area(*hole_it, rect_bbox);
  }
  return area;
}

// ---------------------------------------------------------------------
// Supercover Line and Edge Rasterization
// ---------------------------------------------------------------------

// Compute a supercover line (returns grid cells traversed by the
// line) based on Amanatides–Woo's line algorithm. Perfectly accurate
static GridCoordinates compute_supercover_line_amanatides_woo(
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
  double t_delta_x = almost_equal(delta_x, 0.0)
                       ? std::numeric_limits<double>::infinity()
                       : 1.0 / std::abs(delta_x);
  double t_delta_y = almost_equal(delta_y, 0.0)
                       ? std::numeric_limits<double>::infinity()
                       : 1.0 / std::abs(delta_y);

  // The parametric distance along the ray until we reach the first vertical
  // (or horizontal) grid line
  double t_max_x, t_max_y;
  if (!almost_equal(delta_x, 0.0)) {
    if (step_x > 0)
      t_max_x = ((cell_x + 1) - x0) / delta_x;
    else
      t_max_x = (x0 - cell_x) / -delta_x;
  } else {
    t_max_x = std::numeric_limits<double>::infinity();
  }

  if (!almost_equal(delta_y, 0.0)) {
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

static GridCoordinatesWithEdgeEndpoints rasterize_polygon_edges(
  const Polygon &polygon)
{
  GridCoordinatesWithEdgeEndpoints cells_with_edges;
  if (polygon.size() < 2) {
    return cells_with_edges;
  }

  auto next_pt_ind = [&polygon](unsigned int i) -> unsigned int {
    return (i + 1) % static_cast<unsigned int>(polygon.size());
  };
  for (unsigned int i = 0; i < polygon.size(); ++i) {
    const unsigned int j = next_pt_ind(i);
    const Point p0 = polygon[i];
    const Point p1 = polygon[j];
    GridCoordinates line_cells =
      compute_supercover_line_amanatides_woo(p0.x(), p0.y(), p1.x(), p1.y());
    for (const auto &cell : line_cells) {
      if (cells_with_edges.empty()) {
        cells_with_edges.push_back({cell, {i, i}});
      } else {
        auto &[prev_cell, edge] = cells_with_edges.back();
        if (prev_cell != cell) {
          cells_with_edges.push_back({cell, {i, i}});
        } else {
          edge.second = i;
        }
      }
    }
  }
  // Make last point outside of the cell
  for (auto &[cell, edge] : cells_with_edges) {
    edge.second = next_pt_ind(edge.second);
  }
  // Make sure circular property of the polygon is preserved
  if (
    cells_with_edges.size() > 1 and
    cells_with_edges.front().first == cells_with_edges.back().first) {
    auto &[first_cell, first_edge] = cells_with_edges.front();
    auto &[last_cell, last_edge] = cells_with_edges.back();
    if (first_edge.first == last_edge.second) {
      // Connect the last edge to the first edge and update the first edge
      first_edge.first = last_edge.first;

      // Remove the last edge
      cells_with_edges.pop_back();
    }
  }
  return cells_with_edges;
}

// ---------------------------------------------------------------------
// Point-in-Polygon Check for Polygon_with_holes
// ---------------------------------------------------------------------

static bool is_point_inside_polygon(const Polygon &polygon, const Point &pt)
{
  return polygon.bounded_side(pt) == CGAL::ON_BOUNDED_SIDE;
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
  unsigned int entering_first_edge_idx;
  unsigned int exited_last_edge_idx;
};

struct GeoDivAreaInfo {
  unsigned int gd_id;
  double area;
};

struct Interval {
  unsigned int enter;
  unsigned int exit;
};

// ---------------------------------------------------------------------
// Helper Functions for fill_with_density_clip
// ---------------------------------------------------------------------

// Runs the supercover line algorithm on the edges of the Map
// * For each cell, stores which pwhs' edges are present in the cell
// * Whether the cell is an edge cell
// * Stores metadata for each pwh
static void process_geo_divisions_edge_info(
  boost::multi_array<std::vector<PolygonInfo>, 2> &edge_cell_polyinfo,
  std::vector<PolygonInfo> &all_pwh_info,
  InsetState &inset_state)
{
#pragma omp parallel for collapse(2) schedule(dynamic)
  for (unsigned int gd_id = 0; gd_id < inset_state.geo_divs().size();
       ++gd_id) {
    const GeoDiv &gd = inset_state.geo_divs()[gd_id];
    for (unsigned int pwh_id = 0; pwh_id < gd.n_polygons_with_holes();
         ++pwh_id) {
      const Polygon_with_holes &pwh = gd.polygons_with_holes()[pwh_id];
      const GridCoordinatesWithEdgeEndpoints outer_cells =
        rasterize_polygon_edges(pwh.outer_boundary());
      for (const auto &[cell, poly_idx_pair] : outer_cells) {
        const auto &[x, y] = cell;
        auto &[entering_first_edge_idx, exited_last_edge_idx] = poly_idx_pair;
        const PolygonInfo outer_info{
          gd_id,
          static_cast<unsigned int>(all_pwh_info.size()),
          pwh_id,
          false,
          0,
          entering_first_edge_idx,
          exited_last_edge_idx};
        edge_cell_polyinfo[x][y].push_back(outer_info);
      }

      for (unsigned int hole_id = 0; hole_id < pwh.number_of_holes();
           ++hole_id) {
        const Polygon &hole = pwh.holes()[hole_id];
        const GridCoordinatesWithEdgeEndpoints hole_cells =
          rasterize_polygon_edges(hole);
        for (const auto &[cell, poly_idx_pair] : hole_cells) {
          const auto &[x, y] = cell;
          auto &[entering_first_edge_idx, exited_last_edge_idx] =
            poly_idx_pair;
          const PolygonInfo hole_info{
            gd_id,
            static_cast<unsigned int>(all_pwh_info.size()),
            pwh_id,
            true,
            hole_id,
            entering_first_edge_idx,
            exited_last_edge_idx};
          edge_cell_polyinfo[x][y].push_back(hole_info);
        }
      }
      const PolygonInfo pwh_info{
        gd_id,
        static_cast<unsigned int>(all_pwh_info.size()),
        pwh_id,
        false,
        0,
        0,
        0};
      all_pwh_info.push_back(pwh_info);
    }
  }
}

// For all non-edge cells, compute connected components
// and assign a unique ID to each connected component
// This is useful later to classify all cells of a connected
// component to a single pwh or to ocean
static unsigned int compute_connected_components(
  boost::multi_array<int, 2> &comp,
  const boost::multi_array<std::vector<PolygonInfo>, 2> &edge_cell_polyinfo,
  InsetState &inset_state)
{
  const unsigned int lx = inset_state.lx();
  const unsigned int ly = inset_state.ly();

  // Precompute a boolean mask: true if cell (x,y) is an edge cell
  std::vector<bool> is_edge(lx * ly, false);
  for (unsigned int x = 0; x < lx; ++x) {
    for (unsigned int y = 0; y < ly; ++y) {

      // Flatten index: index = x * ly + y
      is_edge[static_cast<unsigned int>(x * ly + y)] =
        (!edge_cell_polyinfo[x][y].empty());
    }
  }

  unsigned int n_components = 0;
  const int dx[4] = {0, 0, 1, -1};
  const int dy[4] = {1, -1, 0, 0};

  // Flood fill
  for (unsigned int x = 0; x < lx; ++x) {
    for (unsigned int y = 0; y < ly; ++y) {
      if (
        !is_edge[static_cast<unsigned int>(x * ly + y)] && comp[x][y] == -1) {
        std::queue<Coordinate> q;
        q.push({x, y});
        comp[x][y] = static_cast<int>(n_components);
        while (!q.empty()) {
          auto [cx, cy] = q.front();
          q.pop();
          for (int d = 0; d < 4; ++d) {
            const int tx = static_cast<int>(cx) + dx[d];
            const int ty = static_cast<int>(cy) + dy[d];

            if (
              tx < 0 || tx >= static_cast<int>(lx) || ty < 0 ||
              ty >= static_cast<int>(ly))
              continue;

            const unsigned int nx = static_cast<unsigned int>(tx);
            const unsigned int ny = static_cast<unsigned int>(ty);
            if (!is_edge[nx * ly + ny] && comp[nx][ny] == -1) {
              comp[nx][ny] = static_cast<int>(n_components);
              q.push({nx, ny});
            }
          }
        }
        ++n_components;
      }
    }
  }
  return n_components;
}

static Polygon extract_subpath(
  const Polygon &poly,
  unsigned int start_idx,
  unsigned int end_idx)
{
  Polygon subpath;
  const size_t n = poly.size();
  if (n == 0) {
    return subpath;
  }
  if (end_idx > start_idx) {
    for (size_t i = start_idx; i <= end_idx; ++i)
      subpath.push_back(poly[i]);
  } else {
    // Wrap-around: take from start_idx to end, then from 0 to end_idx
    for (size_t i = start_idx; i < n; ++i)
      subpath.push_back(poly[i]);
    for (size_t i = 0; i <= end_idx; ++i)
      subpath.push_back(poly[i]);
  }
  return subpath;
}

// --- Clockwise ordering helpers ---
//
// We define a clockwise ordering of the cell boundary for a 1x1 cell with
// bottom-left (cx,cy): Edge 0 (right): x == cx+1, order from bottom to top.
// Edge 1 (bottom): y == cy, order from right to left.
// Edge 2 (left): x == cx, order from top to bottom.
// Edge 3 (top): y == cy+1, order from left to right.
inline int get_edge_index_cw(
  const Point &p,
  double cx,
  double cy,
  double eps = 1e-7)
{
  if (std::fabs(p.x() - (cx + 1)) < eps)
    return 0;  // right edge
  else if (std::fabs(p.y() - cy) < eps)
    return 1;  // bottom edge
  else if (std::fabs(p.x() - cx) < eps)
    return 2;  // left edge
  else if (std::fabs(p.y() - (cy + 1)) < eps)
    return 3;  // top edge
  else
    throw std::runtime_error("Point is not on the cell boundary.");
}

// Returns the corner for each edge (the extreme point along that edge in the
// direction of traversal). In our clockwise order:
// - On the right edge (0), we traverse upward; the final point is the
// top–right corner.
// - On the bottom edge (1), we traverse from right to left; the first point is
// the right–bottom corner.
// - On the left edge (2), we traverse downward; the final point is the
// bottom–left corner.
// - On the top edge (3), we traverse from left to right; the first point is
// the left–top corner.
inline Point get_corner_cw(int edge, double cx, double cy)
{
  switch (edge) {
  case 0:
    return Point(cx + 1, cy);  // right edge: top–right corner
  case 1:
    return Point(cx, cy);  // bottom edge: right–bottom corner
  case 2:
    return Point(cx, cy + 1);  // left edge: bottom–left corner
  case 3:
    return Point(cx + 1, cy + 1);  // top edge: left–top corner
  default:
    throw std::runtime_error("Invalid edge index.");
  }
}

// Checks the order of two points on the same edge in clockwise order.
// For each edge, the order is defined as follows:
// - Right edge (0): increasing y (from bottom to top).
// - Bottom edge (1): decreasing x (from right to left).
// - Left edge (2): decreasing y (from top to bottom).
// - Top edge (3): increasing x (from left to right).
inline bool on_same_edge_in_order_cw(
  const Point &a,
  const Point &b,
  int edge,
  double eps = 1e-9)
{
  switch (edge) {
  case 0:
    return a.y() > b.y() - eps;  // right edge: from lower to higher y
  case 1:
    return a.x() > b.x() + eps;  // bottom edge: from higher to lower x
  case 2:
    return a.y() < b.y() + eps;  // left edge: from higher to lower y
  case 3:
    return a.x() < b.x() - eps;  // top edge: from lower to higher x
  default:
    throw std::runtime_error("Invalid edge index.");
  }
}

// --- Building the connecting path in clockwise order ---
//
// p1 is the first vertex of the clipped path, p2 is the last vertex.
// Both are assumed to lie exactly on the boundary of a 1x1 cell with
// bottom-left at (cx,cy). We want to “wrap” p2 around the cell boundary in a
// clockwise manner until we reach p1.
static Polygon get_connecting_path_cw(
  const Point &p1,
  const Point &p2,
  double cx,
  double cy)
{
  Polygon connecting_path;
  connecting_path.push_back(p2);

  int edge_p1 = get_edge_index_cw(p1, cx, cy);
  int edge_p2 = get_edge_index_cw(p2, cx, cy);

  // If both endpoints lie on the same edge:
  if (edge_p1 == edge_p2) {
    // On that edge, check if p1 comes after p2 (in our clockwise order)
    if (on_same_edge_in_order_cw(p2, p1, edge_p1)) {
      return connecting_path;
    }
  }

  // If endpoints are on different edges, traverse edges in clockwise order.
  // Our clockwise order is cyclic: 0 -> 1 -> 2 -> 3 -> 0.
  int current_edge = edge_p2;

  // First, complete p2's edge by adding its extreme corner
  Point corner = get_corner_cw(current_edge, cx, cy);
  if (!(corner == p2))
    connecting_path.push_back(corner);

  // Move to the next edge in clockwise order
  current_edge = (current_edge + 1) % 4;

  // Traverse until we reach the edge of p1
  while (current_edge != edge_p1) {
    connecting_path.push_back(get_corner_cw(current_edge, cx, cy));
    current_edge = (current_edge + 1) % 4;
  }
  return connecting_path;
}

// Returns true if point p lies on any edge of the bounding box bbox.
// That is, if p lies on either the top, bottom, left, or right edge.
inline bool is_point_on_rect_edge(
  const Point &p,
  const Bbox &bbox,
  double eps = 1e-9)
{
  bool onBottom = (std::fabs(p.y() - bbox.ymin()) < eps) &&
                  (p.x() >= bbox.xmin() - eps && p.x() <= bbox.xmax() + eps);
  bool onTop = (std::fabs(p.y() - bbox.ymax()) < eps) &&
               (p.x() >= bbox.xmin() - eps && p.x() <= bbox.xmax() + eps);
  bool onLeft = (std::fabs(p.x() - bbox.xmin()) < eps) &&
                (p.y() >= bbox.ymin() - eps && p.y() <= bbox.ymax() + eps);
  bool onRight = (std::fabs(p.x() - bbox.xmax()) < eps) &&
                 (p.y() >= bbox.ymin() - eps && p.y() <= bbox.ymax() + eps);
  return onBottom || onTop || onLeft || onRight;
}

// Checks if the clipped_path is fully contained within the bounding box
inline bool is_path_contained_within_cell(
  const Polygon &clipped_path,
  const Bbox &bbox)
{
  const Point &first = clipped_path[0];
  const Point &last = clipped_path[clipped_path.size() - 1];

  return is_point_on_rect_edge(first, bbox) &&
         is_point_on_rect_edge(last, bbox);
}

// If inside, for each connected component (cc), maps it to corresponding pwh
// For an unmapped cell of a cc, if the neighboring cell is a edge cell,
// then we find the pwh of that edge, and run point-in-polygon check
// to see if the cc is inside the pwh. If so, we map the cc to that pwh
static void map_components_to_pwh(
  std::vector<int> &comp_id_to_pwh_tot_id,
  const boost::multi_array<std::vector<PolygonInfo>, 2> &edge_cell_polyinfo,
  const boost::multi_array<int, 2> &comp,
  const std::vector<PolygonInfo> &all_pwh_info,
  InsetState &inset_state)
{
  const unsigned int lx = inset_state.lx();
  const unsigned int ly = inset_state.ly();

  std::vector<bool> is_edge(lx * ly, false);
  for (unsigned int x = 0; x < lx; ++x) {
    for (unsigned int y = 0; y < ly; ++y) {

      // Flatten index: index = x * ly + y
      is_edge[static_cast<unsigned int>(x * ly + y)] =
        (!edge_cell_polyinfo[x][y].empty());
    }
  }

  auto is_comp_used = [&comp_id_to_pwh_tot_id](int comp_id) {
    return (comp_id != -1) and
           comp_id_to_pwh_tot_id[static_cast<unsigned int>(comp_id)] != -1;
  };

  // Store the connected components that are outside the pwh when encountered
  // for the first time. Avoids running the expensive point-in-polygon check
  // multiple times
  std::vector<std::set<int>> is_outside(all_pwh_info.size());

  const int dx[4] = {0, 0, 1, -1};
  const int dy[4] = {1, -1, 0, 0};

  for (unsigned int x = 0; x < lx; ++x) {
    for (unsigned int y = 0; y < ly; ++y) {
      const int comp_id = comp[x][y];
      if (is_edge[x * ly + y] || is_comp_used(comp_id)) {
        continue;
      }
      for (int d = 0; d < 4; ++d) {
        const int tx = static_cast<int>(x) + dx[d];
        const int ty = static_cast<int>(y) + dy[d];

        if (
          tx < 0 || tx >= static_cast<int>(lx) || ty < 0 ||
          ty >= static_cast<int>(ly))
          continue;

        const unsigned int nx = static_cast<unsigned int>(tx);
        const unsigned int ny = static_cast<unsigned int>(ty);

        if (is_edge[nx * ly + ny]) {
          auto cell_polyinfos = edge_cell_polyinfo[nx][ny];
          for (const auto &poly_info : cell_polyinfos) {
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
              if (is_point_inside_polygon(hole, center)) {
                is_outside[pwh_tot_id].insert(comp_id);
              } else {
                comp_id_to_pwh_tot_id[static_cast<unsigned int>(comp_id)] =
                  static_cast<int>(pwh_tot_id);
              }
            } else {
              auto &outer_boundary = pwh.outer_boundary();
              if (is_point_inside_polygon(outer_boundary, center)) {
                comp_id_to_pwh_tot_id[static_cast<unsigned int>(comp_id)] =
                  static_cast<int>(pwh_tot_id);
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

// Build the closed polygon by appending the connecting path to the clipped
// path. The closed polygon will start at p1 (first vertex of clipped_path),
// then follow the clipped path, then follow the connecting path (which wraps
// clockwise from the last vertex back to p1).
static Polygon build_closed_polygon(
  const Polygon &clipped_path,
  const std::pair<int, int> &cell_bottom_left)
{
  double cx = static_cast<double>(cell_bottom_left.first);
  double cy = static_cast<double>(cell_bottom_left.second);

  // Use indexing to obtain endpoints
  const Point p1 = clipped_path[0];
  const Point p2 = clipped_path[clipped_path.size() - 1];

  // Compute the connecting path along the cell boundary from p1 to p2
  // (clockwise)
  Polygon connecting_path = get_connecting_path_cw(p1, p2, cx, cy);

  Polygon closed_poly = clipped_path;

  // Insert the connecting path, skipping the first point (duplicate of p2)
  closed_poly.insert(
    closed_poly.end(),
    connecting_path.begin() + 1,
    connecting_path.end());

  // Close the polygon by appending p1
  closed_poly.push_back(p1);
  return closed_poly;
}

static double compute_path_union_area(
  const Polygon &polygon,
  bool is_outer,
  const std::vector<Interval> &intervals,
  const Coordinate &cell_coor)
{
  const double cx = cell_coor.first;
  const double cy = cell_coor.second;
  const Bbox cell_bbox(cx, cy, cx + 1, cy + 1);

  Polygon unified_path;

  for (size_t i = 0; i < intervals.size(); ++i) {
    const auto &iv = intervals[i];
    Polygon subpath = extract_subpath(polygon, iv.enter, iv.exit);
    Polygon clipped_path = clip_path_by_rectangle(subpath, cell_bbox);

    // If the first interval's clipped path is not fully contained,
    // fall back to computing the overlap area
    if (i == 0 && !is_path_contained_within_cell(clipped_path, cell_bbox)) {
      return compute_polygon_rectangle_overlap_area(polygon, cell_bbox);
    }

    if (is_outer) {
      // Outer boundaries are oriented in counter-clockwise. We make sure
      // all boundaries are oriented in clockwise order to simplify
      // implementation
      std::reverse(clipped_path.begin(), clipped_path.end());
    }

    if (i == 0) {
      unified_path = clipped_path;
    } else {
      // Connector from the last point of unified_path to the first point of
      // clipped_path
      const Point last_pt = unified_path[unified_path.size() - 1];
      const Point first_pt = clipped_path[0];
      Polygon connecting_path =
        get_connecting_path_cw(first_pt, last_pt, cx, cy);

      // Insert the connecting path, skipping the first point (which duplicates
      // last_pt)
      if (connecting_path.size() > 1) {
        unified_path.insert(
          unified_path.end(),
          connecting_path.begin() + 1,
          connecting_path.end());
      }

      unified_path.insert(
        unified_path.end(),
        clipped_path.begin(),
        clipped_path.end());
    }
  }

  // Close the polygon using the unified path
  Polygon closed_poly = build_closed_polygon(unified_path, cell_coor);
  double area = std::abs(closed_poly.area());
  return area - std::floor(area);
}

static void compute_area(
  const boost::multi_array<std::vector<PolygonInfo>, 2> &edge_cell_polyinfo,
  const boost::multi_array<int, 2> &comp,
  const std::vector<int> &comp_id_to_pwh_tot_id,
  const std::vector<PolygonInfo> &all_pwh_info,
  FTReal2d &rho,
  InsetState &inset_state)
{
  const unsigned int lx = inset_state.lx();
  const unsigned int ly = inset_state.ly();
  std::vector<bool> is_edge(lx * ly, false);
  for (unsigned int x = 0; x < lx; ++x) {
    for (unsigned int y = 0; y < ly; ++y) {

      // Flatten index: index = x * ly + y
      is_edge[static_cast<unsigned int>(x * ly + y)] =
        (!edge_cell_polyinfo[x][y].empty());
    }
  }

  auto is_inside_polygon = [&](unsigned int x, unsigned int y) {
    auto comp_id = comp[x][y];
    return (comp_id != -1) and
           comp_id_to_pwh_tot_id[static_cast<unsigned int>(comp[x][y])] != -1;
  };

  std::vector<double> gd_target_density;
  for (const auto &gd : inset_state.geo_divs()) {
    gd_target_density.push_back(
      inset_state.target_area_at(gd.id()) / gd.area());
  }
  const double ocean_density = (lx * ly - inset_state.total_target_area()) /
                               (lx * ly - inset_state.total_inset_area());
  const double ocean_area_error = std::abs(
    (lx * ly - inset_state.total_inset_area()) /
      (lx * ly - inset_state.total_target_area()) -
    1.0);

#pragma omp parallel for collapse(2) schedule(dynamic)
  for (unsigned int x = 0; x < lx; ++x) {
    for (unsigned int y = 0; y < ly; ++y) {
      const int comp_id = comp[x][y];
      double num = 0.0;
      double den = 0.0;
      double area_tot = 0.0;
      if (is_inside_polygon(x, y)) {
        const unsigned int pwh_tot_id = static_cast<unsigned int>(
          comp_id_to_pwh_tot_id[static_cast<unsigned int>(comp_id)]);
        const PolygonInfo &poly_info = all_pwh_info[pwh_tot_id];
        const unsigned int gd_id = poly_info.gd_id;
        const std::string &gd_name = inset_state.geo_divs()[gd_id].id();
        const double weight = 1.0 * inset_state.area_error_at(gd_name);
        num += weight * gd_target_density[gd_id];
        den += weight;
        area_tot += 1.0;
      } else if (is_edge[static_cast<unsigned int>(x * ly + y)]) {
        auto cell_poly_info = edge_cell_polyinfo[x][y];
        const size_t n = cell_poly_info.size();

        for (unsigned int i = 0; i < n;) {
          auto &poly_info = cell_poly_info[i];
          const unsigned int pwh_tot_id = poly_info.pwh_tot_id;
          double outer_area = 1.0;
          double hole_taken_area = 0.0;
          if (not poly_info.is_hole) {
            unsigned int j = i;
            while (j + 1 < n and
                   cell_poly_info[j + 1].pwh_tot_id == pwh_tot_id and
                   not cell_poly_info[j + 1].is_hole) {
              ++j;
            }
            const auto &polygon = inset_state.geo_divs()[poly_info.gd_id]
                                    .polygons_with_holes()[poly_info.pwh_id]
                                    .outer_boundary();
            std::vector<Interval> intervals;
            for (unsigned int k = i; k <= j; ++k) {
              const unsigned int enter =
                cell_poly_info[k].entering_first_edge_idx;
              const unsigned int exit = cell_poly_info[k].exited_last_edge_idx;
              intervals.push_back({enter, exit});
            }
            outer_area =
              compute_path_union_area(polygon, true, intervals, {x, y});
            i = j + 1;
          }
          while (i < n and cell_poly_info[i].pwh_tot_id == pwh_tot_id) {
            const unsigned int hole_id = cell_poly_info[i].hole_id;
            const auto &polygon = inset_state.geo_divs()[poly_info.gd_id]
                                    .polygons_with_holes()[poly_info.pwh_id]
                                    .holes()[hole_id];
            unsigned int j = i;
            while (j + 1 < n and
                   cell_poly_info[j + 1].pwh_tot_id == pwh_tot_id and
                   cell_poly_info[j + 1].hole_id == hole_id) {
              ++j;
            }
            std::vector<Interval> intervals;
            for (unsigned int k = i; k <= j; ++k) {
              const unsigned int enter =
                cell_poly_info[k].entering_first_edge_idx;
              const unsigned int exit = cell_poly_info[k].exited_last_edge_idx;
              intervals.push_back({enter, exit});
            }
            hole_taken_area +=
              compute_path_union_area(polygon, false, intervals, {x, y});
            i = j + 1;
          }
          const double area = outer_area - hole_taken_area;
          const unsigned int gd_id = poly_info.gd_id;
          const std::string &gd_name = inset_state.geo_divs()[gd_id].id();
          const double weight = area * inset_state.area_error_at(gd_name);
          num += weight * gd_target_density[gd_id];
          den += weight;
          area_tot += area;
        }
      }

      const double ocean_area = 1.0 - area_tot;
      const double ocean_weight = ocean_area * ocean_area_error;
      num += ocean_weight * ocean_density;
      den += ocean_weight;

      if (den > 0.0) {
        rho(x, y) = num / den;
      } else {
        rho(x, y) = ocean_density;
      }
    }
  }
}

// Test function to compare the computed area with the actual area (brute
// force)
[[maybe_unused]] static void test_areas_densities(
  const FTReal2d &rho,
  InsetState &inset_state)
{
  const unsigned int lx = inset_state.lx();
  const unsigned int ly = inset_state.ly();
  boost::multi_array<double, 2> numer(boost::extents[lx][ly]);
  boost::multi_array<double, 2> denom(boost::extents[lx][ly]);
  boost::multi_array<double, 2> area_found(boost::extents[lx][ly]);

  std::fill_n(numer.data(), numer.num_elements(), 0.0);
  std::fill_n(denom.data(), denom.num_elements(), 0.0);
  std::fill_n(area_found.data(), area_found.num_elements(), 0.0);

  std::vector<double> gd_target_density;
  for (const auto &gd : inset_state.geo_divs()) {
    gd_target_density.push_back(
      inset_state.target_area_at(gd.id()) / gd.area());
  }

  boost::multi_array<double, 2> area_actual(boost::extents[lx][ly]);
  std::fill_n(area_actual.data(), area_actual.num_elements(), 0.0);
  for (unsigned int gd_id = 0; gd_id < inset_state.geo_divs().size();
       ++gd_id) {
    const GeoDiv &gd = inset_state.geo_divs()[gd_id];
    for (const auto &pwh : gd.polygons_with_holes()) {
      auto bbox = pwh.outer_boundary().bbox();
      for (unsigned int i = static_cast<unsigned int>(
             std::max(0, static_cast<int>(bbox.xmin())));
           i < std::min(lx, static_cast<unsigned int>(bbox.xmax()) + 1);
           ++i) {
        for (unsigned int j = static_cast<unsigned int>(
               std::max(0, static_cast<int>(bbox.ymin())));
             j < std::min(ly, static_cast<unsigned int>(bbox.ymax()) + 1);
             ++j) {

          // Create a 1x1 cell
          Polygon cell;
          cell.push_back(Point(i, j));
          cell.push_back(Point(i + 1, j));
          cell.push_back(Point(i + 1, j + 1));
          cell.push_back(Point(i, j + 1));

          const double intersect_area_pwh =
            compute_pwh_rectangle_overlap_area(pwh, cell.bbox());
          area_actual[i][j] += intersect_area_pwh;
          const double weight =
            intersect_area_pwh * inset_state.area_error_at(gd.id());
          numer[i][j] += weight * gd_target_density[gd_id];
          denom[i][j] += weight;
        }
      }
    }
  }

  const double ocean_density = (lx * ly - inset_state.total_target_area()) /
                               (lx * ly - inset_state.total_inset_area());
  const double ocean_area_error = std::abs(
    (lx * ly - inset_state.total_inset_area()) /
      (lx * ly - inset_state.total_target_area()) -
    1.0);

  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      double weight = (1.0 - area_actual[i][j]) * ocean_area_error;
      numer[i][j] += weight * ocean_density;
      denom[i][j] += weight;
    }
  }

  FTReal2d rho_actual;
  rho_actual.allocate(lx, ly);

  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      if (denom[i][j] > 0)
        rho_actual(i, j) = numer[i][j] / denom[i][j];
      else
        rho_actual(i, j) = ocean_density;
    }
  }

  double absolute_error_density = 0.0;
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      absolute_error_density += std::abs(rho_actual(i, j) - rho(i, j));
    }
  }
  const unsigned int num_cells = lx * ly;
  const double mean_absolute_error_density =
    absolute_error_density / num_cells;
  std::cout << "Num Integration: " << inset_state.n_finished_integrations()
            << std::endl;
  std::cout << "Mean Absolute Error Density: " << mean_absolute_error_density
            << std::endl;
}

void InsetState::create_contiguity_graph()
{
  // Create a copy of `inset_state` to avoid rescaling the original
  InsetState is_copy = *this;

  is_copy.rescale_map();

  boost::multi_array<std::vector<PolygonInfo>, 2> edge_cell_polyinfo(
    boost::extents[is_copy.lx_][is_copy.ly_]);
  std::vector<PolygonInfo> all_pwh_info;

  process_geo_divisions_edge_info(edge_cell_polyinfo, all_pwh_info, is_copy);

  for (unsigned int x = 0; x < is_copy.lx(); ++x) {
    for (unsigned int y = 0; y < is_copy.ly(); ++y) {
      std::vector<PolygonInfo> polygons_at_xy = edge_cell_polyinfo[x][y];

      for (size_t i = 0; i < polygons_at_xy.size(); ++i) {
        GeoDiv &gd_1 = geo_divs_[polygons_at_xy[i].gd_id];
        for (size_t j = i + 1; j < polygons_at_xy.size(); ++j) {
          GeoDiv &gd_2 = geo_divs_[polygons_at_xy[j].gd_id];

          // Update adjacency for both GeoDivs
          // GeoDiv.adjacent_to() already checks whether the GeoDivs are
          // different
          gd_1.adjacent_to(gd_2.id());
          gd_2.adjacent_to(gd_1.id());
        }
      }
    }
  }
}

void InsetState::fill_with_density_clip()
{
  std::cerr << "Filling density using clipping method" << std::endl;

  timer.start("Fill with Density (Clipping Method)");

  // Step 1: Detect edges and store edge information
  boost::multi_array<std::vector<PolygonInfo>, 2> edge_cell_polyinfo(
    boost::extents[lx_][ly_]);
  std::vector<PolygonInfo> all_pwh_info;

  process_geo_divisions_edge_info(edge_cell_polyinfo, all_pwh_info, *this);

  // Step 2: Compute connected components
  boost::multi_array<int, 2> comp(boost::extents[lx_][ly_]);
  std::fill_n(comp.data(), comp.num_elements(), -1);
  const unsigned int n_components =
    compute_connected_components(comp, edge_cell_polyinfo, *this);

  // Step 3: Map components to polygon with holes IDs
  std::vector<int> comp_id_to_pwh_tot_id(n_components, -1);
  map_components_to_pwh(
    comp_id_to_pwh_tot_id,
    edge_cell_polyinfo,
    comp,
    all_pwh_info,
    *this);

  compute_area(
    edge_cell_polyinfo,
    comp,
    comp_id_to_pwh_tot_id,
    all_pwh_info,
    rho_init_,
    *this);

  // test_areas_densities(rho_init_, *this);

  // Step 6: Run FFTW to compute rho_ft_
  auto rho_begin = rho_init_.as_1d_array();
  auto rho_end = rho_begin + lx_ * ly_;
  auto [min_iter, max_iter] = std::minmax_element(rho_begin, rho_end);
  dens_min_ = *min_iter;
  dens_max_ = *max_iter;

  execute_fftw_fwd_plan();
  timer.stop("Fill with Density (Clipping Method)");
}
