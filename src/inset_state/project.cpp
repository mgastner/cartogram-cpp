#include "inset_state.hpp"
#include "interpolate_bilinearly.hpp"
#include "round_point.hpp"
#include "segments_that_meet.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_set>

bool InsetState::project()
{

  if (!create_delaunay_t())  // Triangle has flipped during triangulation
    return false;

  if (!args_.disable_simplification_densification) {
    densify_geo_divs_using_delaunay_t();
  }

  // Plot if requested
  if (args_.plot_quadtree) {
    write_delaunay_triangles(
      file_prefix_ + "c_updated_delaunay_t_after_flatten",
      false);
  }

  // Project using the updated Delaunay triangulation and plot
  project_with_delaunay_t(args_.redirect_exports_to_stdout);

  if (args_.plot_quadtree) {
    write_delaunay_triangles(
      file_prefix_ + "d_projected_with_updated_delaunay_t",
      true);
  }

  if (!args_.disable_simplification_densification) {

    simplify(args_.target_points_per_inset);
  }
  if (args_.plot_intersections) {
    write_intersections_image();
  }

  return true;
}

template <typename Triangle>
static Point interpolate_point_with_barycentric_coordinates(
  const Point &p,
  const Triangle &triangle,
  const ProjectionData &proj)
{
  // Get the three vertices
  const Point v1 = triangle.vertices[0];
  const Point v2 = triangle.vertices[1];
  const Point v3 = triangle.vertices[2];

  // Calculate barycentric coordinates
  auto bary =
    CGAL::Barycentric_coordinates ::triangle_coordinates_in_tuple_2<Point>(
      v1,
      v2,
      v3,
      p);

  // Get the barycentric coordinates
  const double bary_x = std::get<0>(bary);
  const double bary_y = std::get<1>(bary);
  const double bary_z = std::get<2>(bary);

  auto to_uint = [](const double val) {
    return static_cast<uint32_t>(val + 0.5);
  };

  // Get projected vertices
  const Point v1_proj = proj.get(to_uint(v1.x()), to_uint(v1.y()));
  const Point v2_proj = proj.get(to_uint(v2.x()), to_uint(v2.y()));
  const Point v3_proj = proj.get(to_uint(v3.x()), to_uint(v3.y()));

  // Calculate projected point of p
  return {
    bary_x * v1_proj.x() + bary_y * v2_proj.x() + bary_z * v3_proj.x(),
    bary_x * v1_proj.y() + bary_y * v2_proj.y() + bary_z * v3_proj.y()};
}

void InsetState::project_with_delaunay_t(bool output_to_stdout)
{
  timer.start("Project");
  auto lambda_bary = [&](Point p1) {
    return interpolate_point_with_barycentric_coordinates(
      p1,
      *triang_.locate(p1),
      proj_data_);
  };
  transform_points(lambda_bary);

  if (output_to_stdout) {
    transform_points(lambda_bary, true);
  }

  // Every ring passes the test of is_simple() unless the repair returns
  // false; is_simple() then reports the first ring that fails
  if (
    !args_.disable_simplification_densification &&
    !remove_vertices_where_rings_cross())
    is_simple(__func__);
  timer.stop("Project");
}

// Distance from p to the segment from a to b
static double distance_to_segment(
  const Point &p,
  const Point &a,
  const Point &b)
{
  const double dx = b.x() - a.x();
  const double dy = b.y() - a.y();
  const double l2 = dx * dx + dy * dy;
  double t = 0.0;
  if (l2 > 0.0) {
    t = ((p.x() - a.x()) * dx + (p.y() - a.y()) * dy) / l2;
    t = std::clamp(t, 0.0, 1.0);
  }
  return std::hypot(p.x() - (a.x() + t * dx), p.y() - (a.y() + t * dy));
}

// The point to remove where two segments of a ring meet: for neighbors that
// fold back, their common end; for segments that touch only at an end that
// occurs twice in the ring, that end; otherwise, of the four ends, the one
// nearest to the other segment
static Point point_to_remove(const Polygon &ring, const SegmentsThatMeet &meet)
{
  const std::size_t n = ring.size();
  auto next = [n](std::size_t k) {
    return (k + 1 == n) ? 0 : k + 1;
  };
  if (meet.fold_back)
    return (next(meet.i) == meet.j) ? ring[meet.j] : ring[meet.i];
  if (meet.touch_at_common_end)
    return ring[meet.u];
  const std::size_t ends[4] = {meet.i, next(meet.i), meet.j, next(meet.j)};
  const std::size_t others[4] = {meet.j, meet.j, meet.i, meet.i};
  std::size_t nearest = 0;
  double nearest_distance = dbl_inf;
  for (std::size_t k = 0; k < 4; ++k) {
    const double d = distance_to_segment(
      ring[ends[k]],
      ring[others[k]],
      ring[next(others[k])]);
    if (d < nearest_distance) {
      nearest = k;
      nearest_distance = d;
    }
  }
  return ring[ends[nearest]];
}

// The projection maps each point with the affine map of the triangle that
// contains it and rounds the result, which can move the point by about a
// unit in the last place of its coordinates. Where a vertex lies closer than
// that to another segment of its ring, for example where earlier
// integrations squeezed a part of a polygon to below the resolution of
// doubles, rounding can put the vertex across that segment, so that the ring
// crosses itself, although no triangle of the projection is inverted. Every
// ring is therefore checked with the simplicity test of the Polygon type,
// which is also the test of InsetState::is_simple(). Where a ring fails it,
// the segments that meet are found with exact predicates, and for each pair,
// one point is removed: the end nearest to the other segment, the common end
// of neighbors that fold back, or the end where two segments touch. A point
// is removed from every ring that contains it, so that shared boundaries stay
// identical, unless fewer than three points would remain. The rings that lost
// points are checked again, with the test of the Polygon type and, if they
// passed before, with exact predicates, until no ring changes. Returns
// whether every ring passes the test of the Polygon type.
bool InsetState::remove_vertices_where_rings_cross()
{
  std::vector<Polygon *> rings;
  for (auto &gd : geo_divs_) {
    for (auto &pwh : gd.ref_to_polygons_with_holes()) {
      rings.push_back(&pwh.outer_boundary());
      for (auto &h : pwh.holes()) {
        rings.push_back(&h);
      }
    }
  }

  std::vector<bool> fails(rings.size(), false);
  std::vector<bool> changed(rings.size(), false);

  // Whether each changed ring was simple with exact predicates before it
  // changed
  std::vector<bool> exactly_simple_before(rings.size(), false);
  std::vector<std::size_t> to_check(rings.size());
  std::iota(to_check.begin(), to_check.end(), std::size_t{0});
  std::unordered_set<Point> removed;
  while (!to_check.empty()) {
    std::unordered_set<Point> drop;
    for (const std::size_t r : to_check) {
      const Polygon &ring = *rings[r];
      fails[r] = !ring.is_simple();
      const bool repair =
        fails[r] || (changed[r] && exactly_simple_before[r] &&
                     !is_simple_with_exact_predicates(ring));
      if (!repair)
        continue;
      for (const SegmentsThatMeet &meet : segments_that_meet(ring)) {
        drop.insert(point_to_remove(ring, meet));
      }
    }

    if (drop.empty())
      break;

    to_check.clear();
    for (std::size_t r = 0; r < rings.size(); ++r) {
      Polygon &ring = *rings[r];
      Polygon kept;
      for (const Point &p : ring) {
        if (drop.count(p) == 0)
          kept.push_back(p);
      }
      if (kept.size() == ring.size() || kept.size() < 3)
        continue;
      if (!changed[r]) {
        exactly_simple_before[r] = is_simple_with_exact_predicates(ring);
        changed[r] = true;
      }
      for (const Point &p : ring) {
        if (drop.count(p) > 0)
          removed.insert(p);
      }
      ring = std::move(kept);
      to_check.push_back(r);
    }
  }
  if (!removed.empty()) {
    std::cerr << "Removed " << removed.size()
              << " point(s) where the projection made a ring cross itself"
              << std::endl;
  }
  return std::none_of(fails.begin(), fails.end(), [](bool f) {
    return f;
  });
}

// Apply projection to all points in set
void InsetState::project_point_set(std::unordered_set<Point> &unprojected)
{
  auto lambda_bary = [&](Point p1) {
    return interpolate_point_with_barycentric_coordinates(
      p1,
      *triang_.locate(p1),
      proj_data_);
  };
  std::unordered_set<Point> projected;
  for (const Point &pt : unprojected) {
    Point pp = lambda_bary(pt);
    projected.insert(pp);
  }
  unprojected = std::move(projected);
}
