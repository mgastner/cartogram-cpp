#include "inset_state.hpp"
#include "round_point.hpp"
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <variant>

// For printing a vector (debugging purposes)
template <typename A>
std::ostream &operator<<(std::ostream &cout, std::vector<A> const &v)
{
  cout << "[";
  for (int i = 0; i < int(v.size()); i++) {
    if (i)
      cout << ", ";
    cout << v[i];
  }
  return cout << "]";
}

// A point location at (-1, -1) is a sign that a point is not on the
// [0, lx]-by-[0, ly] grid used for calculating the density to be equalized

// TODO: Is there a way to avoid #define? The apparent alternative
// `constexpr out_of_range(-1.0, -1.0);` does not work because "the
// type 'const Point' ... is not literal".
#define OUT_OF_RANGE Point(-1.0, -1.0)

// This function has the following parameters:
// - a segment defined by points a and b.
// - a line defined by coef_x, coef_y, and coef_const with the formula
//      coef_x * x + coef_y * y + coef_const = 0.
// The function returns the unique intersection point between them. If this
// intersection point does not exist, the function returns the point called
// OUT_OF_RANGE, which is always outside any grid cell.
static Point calc_intersection(
  const Point &a,
  const Point &b,
  const double coef_x,
  const double coef_y,
  const double coef_const)
{
  const auto result =
    CGAL::intersection(Line(coef_x, coef_y, coef_const), Segment(a, b));
  if (result) {

    // The result of CGAL::intersection can either be a segment, a point, or
    // null. We only want point intersections. Where there is no point
    // intersection, we get a null pointer.
    if (auto p = std::get_if<Point>(&*result))
      return *p;
  }
  return OUT_OF_RANGE;
}

// This function takes points `a` and `b` as arguments, which define a
// segment. The function also takes the following arguments, which define the
// types of diagonals for which we will calculate intersections:
// - slope: the slope of every diagonal.
// - base_intercept: the intercept that is the closest to 0 of a diagonal
//   on the grid. This value is either 0, 0.25, or 0.5.
// - step: what we need to add to each diagonal's intercept to obtain the
//   next diagonal.
using PointLess = bool (*)(const Point &, const Point &);

static void add_diag_inter(
  std::set<Point, PointLess> *intersections,
  const Point &a,
  const Point &b,
  double slope,
  double base_intercept,
  double step,
  const unsigned int lx,
  const unsigned int ly)
{
  double intercept_start =
    floor(std::min(a.y() - slope * a.x(), b.y() - slope * b.x())) +
    base_intercept;
  double intercept_end =
    std::max(a.y() - slope * a.x(), b.y() - slope * b.x());
  for (double d = intercept_start; d <= intercept_end; d += step) {

    // We consider four types of diagonals:
    // - 'normal': y = x + d,
    // - 'antinormal': y = -x + d,
    // - 'steep': y = 2x + d + base_intercept
    //      (where base_intercept = 0.5),
    // - 'antisteep': y = -2x + d + base_intercept
    //      (where base_intercept = 0.5),
    // - 'gentle': y = 0.5x + d + base_intercept,
    //      (where base_intercept = 0.25),
    // - 'antigentle': y = -0.5x + d + base_intercept,
    //      (where base_intercept = 0.25),
    // where d is the double increasing by `step` after each iteration of the
    // loop.
    // Steep and antisteep diagonals appear in grid cells near x = 0
    // and x = lx. Gentle and antigentle diagonals appear in grids near
    // y = 0 and y = ly.
    Point inter = calc_intersection(a, b, slope, -1.0, d);
    bool on_left_or_right_edge = inter.x() < 0.5 || inter.x() > (lx - 0.5);
    bool on_top_or_bottom_edge = inter.y() < 0.5 || inter.y() > (ly - 0.5);
    const double abs_slope = std::fabs(slope);
    if (
      inter != OUT_OF_RANGE &&
      ((almost_equal(abs_slope, 2.0) && on_left_or_right_edge) ||
       (almost_equal(abs_slope, 0.5) && on_top_or_bottom_edge) ||
       (almost_equal(abs_slope, 1) &&
        (on_left_or_right_edge == on_top_or_bottom_edge)))) {
      (*intersections).insert(inter);
    }
  }
}

// TODO: If a or b are themselves intersection points (e.g., if pt1.x is an
// integer plus 0.5), it appears to be included in the returned intersections.
// Would this property cause the point to be included twice in the line
// segment (once when the end point is the argument pt1 and a second time when
// the same end point is the argument pt2)?

// This function takes two points (called pt1 and pt2) and returns all
// horizontal and vertical intersections of the line segment between pt1 and
// pt2 with a grid whose grid lines are placed one unit apart. The
// function also returns all intersections with the diagonals of these
// grid cells. The function assumes that grid cells start at
// (0.5, 0.5).
static std::vector<Point> densification_points(
  const Point &pt1,
  const Point &pt2,
  const unsigned int lx,
  const unsigned int ly)
{
  // If the input points are identical, return them without calculating
  // intersections
  if (almost_equal(pt1, pt2)) {
    std::vector<Point> points;
    points.push_back(pt1);
    points.push_back(pt2);
    return points;
  }

  // Ordered set for storing intersections before removing duplicates
  std::set<Point, PointLess> temp_intersections(less_than);

  // Store the leftmost point of p1 and pt2 as `a`. If both points have the
  // same x-coordinate, then store the lower point as `a`. The other point is
  // stored as `b`. The segments (a, b) and (b, a) describe the same segment.
  // However, if we flip the order of a and b, the resulting intersections are
  // not necessarily the same because of floating point errors.

  // TODO: IN THE COMMENT ABOVE, DOES THE REMARK ABOUT THE FLOATING-POINT
  // ERRORS STILL APPLY AFTER SWITCHING TO SIMPLE CARTESIAN COORDINATES?
  Point a, b;
  if (
    less_than(pt2.x(), pt1.x()) ||
    (almost_equal(pt1.x(), pt2.x()) && less_than(pt2.y(), pt1.y()))) {
    a = pt2;
    b = pt1;
  } else {
    a = pt1;
    b = pt2;
  }
  temp_intersections.insert(a);
  temp_intersections.insert(b);

  // Get vertical intersections
  double x_start = floor(a.x() + 0.5) + 0.5;
  double x_end = b.x();
  for (double x = x_start; x <= x_end; x += almost_equal(x, 0.0) ? 0.5 : 1.0) {
    Point inter = calc_intersection(a, b, 1.0, 0.0, -x);
    if (inter != OUT_OF_RANGE) {
      temp_intersections.insert(inter);
    }
  }

  // Get horizontal intersections
  double y_start = floor(std::min(a.y(), b.y()) + 0.5) + 0.5;
  double y_end = std::max(a.y(), b.y());
  for (double y = y_start; y <= y_end; y += almost_equal(y, 0.0) ? 0.5 : 1.0) {
    Point inter = calc_intersection(a, b, 0.0, 1.0, -y);
    if (inter != OUT_OF_RANGE) {
      temp_intersections.insert(inter);
    }
  }

  // Get bottom-left to top-right diagonal intersections
  add_diag_inter(&temp_intersections, a, b, 1.0, 0.0, 1.0, lx, ly);

  // Get top-left to bottom-right diagonal intersections
  add_diag_inter(&temp_intersections, a, b, -1.0, 0.0, 1.0, lx, ly);

  // Add edge diagonals when at least one point is near the edge of the grid
  if (a.x() < 0.5 || b.x() < 0.5 || a.x() > (lx - 0.5) || b.x() > (lx - 0.5)) {

    // Bottom-left to top-right edge diagonals
    add_diag_inter(&temp_intersections, a, b, 2.0, 0.5, 1.0, lx, ly);

    // Top-left to bottom-right edge diagonals
    add_diag_inter(&temp_intersections, a, b, -2.0, 0.5, 1.0, lx, ly);
  }
  if (a.y() < 0.5 || b.y() < 0.5 || a.y() > (ly - 0.5) || b.y() > (ly - 0.5)) {

    // Bottom-left to top-right edge diagonals
    add_diag_inter(&temp_intersections, a, b, 0.5, 0.25, 0.5, lx, ly);

    // Top-left to bottom-right edge diagonals
    add_diag_inter(&temp_intersections, a, b, -0.5, 0.25, 0.5, lx, ly);
  }

  // Create a Point vector from the set
  std::vector<Point> intersections(
    temp_intersections.begin(),
    temp_intersections.end());

  // Reverse if needed
  if (
    less_than(pt2.x(), pt1.x()) ||
    (almost_equal(pt1.x(), pt2.x()) && less_than(pt2.y(), pt1.y()))) {
    std::reverse(intersections.begin(), intersections.end());
  }
  return intersections;
}

void InsetState::densify_geo_divs()
{
  timer.start("Densification (using Grid Diagonals)");
  std::cerr << "Densifying" << std::endl;
  std::vector<GeoDiv> geodivs_dens;
  for (const auto &gd : geo_divs_) {
    GeoDiv gd_dens(gd.id());
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto &outer = pwh.outer_boundary();
      Polygon outer_dens;

      // Iterate over each point in the outer boundary of the polygon
      for (unsigned int i = 0; i < outer.size(); ++i) {

        // The segment defined by points `a` and `b` is to be densified.
        // `b` should be the vertex of the boundary immediately after `a`,
        // unless `a` is the final vertex of the boundary, in which case `b`
        // should be the first vertex.
        const auto a = outer[i];
        const auto b = (i == outer.size() - 1) ? outer[0] : outer[i + 1];

        // Densify the segment
        const std::vector<Point> outer_pts_dens =
          densification_points(a, b, lx_, ly_);

        // Push all points. Omit the last point because it will be included
        // in the next iteration. Otherwise, we would have duplicated points
        // in the polygon.
        for (unsigned int j = 0; j < (outer_pts_dens.size() - 1); ++j) {
          outer_dens.push_back(outer_pts_dens[j]);
        }
      }
      std::vector<Polygon> holes_v_dens;

      // Iterate over each hole
      for (auto const &h : pwh.holes()) {
        Polygon hole_dens;
        for (unsigned int j = 0; j < h.size(); ++j) {

          // `c` and `d` are determined in the same way as `a` and `b` above
          const Point c = h[j];
          const Point d = (j == h.size() - 1) ? h[0] : h[j + 1];
          const std::vector<Point> hole_pts_dens =
            densification_points(c, d, lx_, ly_);
          for (unsigned int i = 0; i < (hole_pts_dens.size() - 1); ++i) {
            hole_dens.push_back(hole_pts_dens[i]);
          }
        }
        holes_v_dens.push_back(hole_dens);
      }
      const Polygon_with_holes pwh_dens(
        outer_dens,
        holes_v_dens.begin(),
        holes_v_dens.end());
      gd_dens.push_back(pwh_dens);
    }
    geodivs_dens.push_back(gd_dens);
  }
  geo_divs_ = std::move(geodivs_dens);

  is_simple(__func__);
  timer.stop("Densification (using Grid Diagonals)");
}

static std::vector<Point> densification_points_with_delaunay_t(
  const Point &pt1,
  const Point &pt2,
  const Delaunay &dt)
{
  std::vector<Point> dens_points;
  dens_points.reserve(8);

  // keep track of visited segments
  std::unordered_set<Segment> vis_seg;

  // If the input points are identical, return them without calculating
  // intersections
  if (almost_equal(pt1, pt2)) {
    return {pt1, pt2};
  }

  Face_handle f1 = dt.locate(pt1);
  Face_handle f2 = dt.locate(pt2);

  // If they are inside the same triangle then return original points
  if (f1 == f2) {
    return {pt1, pt2};
  }

  Segment segment(pt1, pt2);  // segment to be densified

  Line_face_circulator lfc = dt.line_walk(pt1, pt2);
  Line_face_circulator lfc_begin = lfc;  // we store the begining iterator

  // keeping track of states
  bool f1_first = false, found_a_face = false;
  if (lfc_begin != nullptr) {
    do {
      Face_handle fh = lfc;

      if (found_a_face) {
        if (f1 == fh or f2 == fh) {  // if current face is second face
          break;
        }
      } else if (fh == f1) {  // have not found a matching face yet; we start
                              // with pt1 point
        f1_first = true;
        found_a_face = true;

        // add the pt1 to the list of densification points
        dens_points.push_back(pt1);
      } else if (fh == f2) {  // we start with pt2 point
        found_a_face = true;

        // add the pt2 to the list of densification points
        dens_points.push_back(pt2);
      } else {  // the face does not intersect the segment
        ++lfc;
        continue;
      }

      // create three segments from the triangle
      Segment s1(fh->vertex(0)->point(), fh->vertex(1)->point());
      Segment s2(fh->vertex(1)->point(), fh->vertex(2)->point());
      Segment s3(fh->vertex(2)->point(), fh->vertex(0)->point());

      for (const Segment &tri_seg : {s1, s2, s3}) {
        Segment tri_seg_rev(tri_seg.target(), tri_seg.source());

        // if the segment (undirected) is already visited, continue
        if (
          vis_seg.find(tri_seg) != vis_seg.end() ||
          vis_seg.find(tri_seg_rev) != vis_seg.end()) {
          continue;
        }

        // Update the visited segments
        vis_seg.insert(tri_seg);
        vis_seg.insert(tri_seg_rev);

        Point pt_intersec;
        if (CGAL::do_intersect(segment, tri_seg)) {
          CGAL::Object p = CGAL::intersection(segment, tri_seg);

          // The CGAL::Object p does not need to be a point, but could also be
          // empty or a line. CGAL::assign() only returns `true` if p is a
          // point.
          if (CGAL::assign(pt_intersec, p)) {
            dens_points.push_back(pt_intersec);
          }
        }
      }

      // move the iterator to the next one
      ++lfc;
    } while (lfc != lfc_begin);
  }

  if (dens_points.size() <= 1) {
    return {pt1, pt2};
  }

  if (f1_first) {
    dens_points.push_back(pt2);
  } else {
    dens_points.push_back(pt1);
  }

  // if densification points are in reverse order, reverse them
  if (!almost_equal(dens_points.front(), pt1)) {
    std::reverse(dens_points.begin(), dens_points.end());
  }

  // check validity of densification points: in case the first and last
  // points are not the originally given points, we consider the densificaiton
  // points invalid and return the original points
  if (
    !almost_equal(dens_points.front(), pt1) ||
    !almost_equal(dens_points.back(), pt2)) {
    return {pt1, pt2};
  }

  // keep densified points those are different from the first and last points
  std::vector<Point> dens_points_unique;
  dens_points_unique.reserve(dens_points.size());
  dens_points_unique.push_back(dens_points.front());
  for (unsigned int i = 1; i < dens_points.size() - 1; ++i) {
    if (
      !almost_equal(dens_points[i], dens_points.front()) &&
      !almost_equal(dens_points[i], dens_points.back()) &&
      !almost_equal(dens_points[i], dens_points[i - 1])) {
      dens_points_unique.push_back(dens_points[i]);
    }
  }
  dens_points_unique.push_back(dens_points.back());

  return dens_points_unique;
}

void InsetState::densify_geo_divs_using_delaunay_t()
{
  timer.start("Densification (using Delanuay Triangles)");
  // Empty set containing densified points from previous iteration
  // points_from_densification_.clear();
  // points_before_densification_.clear();

  std::cerr << "Densifying using Delaunay Triangulation" << std::endl;
  size_t n_points_bef_densification = n_points();
  std::cerr << "Points Before Densification: " << n_points_bef_densification
            << std::endl;
  std::vector<GeoDiv> geodivs_dens;
  geodivs_dens.reserve(geo_divs_.size());

  for (const auto &gd : geo_divs_) {
    GeoDiv gd_dens(gd.id());

    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto &outer = pwh.outer_boundary();
      Polygon outer_dens;
      outer_dens.reserve(outer.size() * 5);

      // Iterate over each point in the outer boundary of the polygon
      const std::size_t outer_sz = outer.size();
      for (std::size_t i = 0; i < outer_sz; ++i) {

        // The segment defined by points `a` and `b` is to be densified.
        // `b` should be the vertex of the boundary immediately after `a`,
        // unless `a` is the final vertex of the boundary, in which case `b`
        // should be the first vertex.
        const Point a = outer[i];
        const Point b = (i + 1 == outer_sz) ? outer[0] : outer[i + 1];

        // Densify the segment
        const std::vector<Point> outer_pts_dens =
          densification_points_with_delaunay_t(a, b, proj_qd_.dt);

        // Push all points. Omit the last point because it will be included
        // in the next iteration. Otherwise, we would have duplicated points
        // in the polygon.
        if (outer_pts_dens.size() > 1)
          outer_dens.insert(
            outer_dens.end(),
            outer_pts_dens.begin(),
            outer_pts_dens.end() - 1);
      }

      // Add new points to points_added_via_densification for plotting
      // std::unordered_set<Point> new_in_outer = new_points(outer,
      // outer_dens); points_from_densification_.insert(
      //   new_in_outer.begin(),
      //   new_in_outer.end());
      // points_before_densification_.insert(outer.begin(), outer.end());

      std::vector<Polygon> holes_v_dens;
      holes_v_dens.reserve(pwh.number_of_holes());

      // Iterate over each hole
      for (auto const &h : pwh.holes()) {

        Polygon hole_dens;
        hole_dens.reserve(h.size() * 2);

        const std::size_t h_sz = h.size();
        for (std::size_t j = 0; j < h_sz; ++j) {

          // `c` and `d` are determined in the same way as `a` and `b` above
          const Point c = h[j];
          const Point d = (j + 1 == h_sz) ? h[0] : h[j + 1];

          const std::vector<Point> hole_pts_dens =
            densification_points_with_delaunay_t(c, d, proj_qd_.dt);

          if (hole_pts_dens.size() > 1)
            hole_dens.insert(
              hole_dens.end(),
              hole_pts_dens.begin(),
              hole_pts_dens.end() - 1);
        }

        // Add new points to points_added_via_densification for plotting
        // std::unordered_set<Point> new_in_hole = new_points(h, hole_dens);
        // points_from_densification_.insert(
        //   new_in_hole.begin(),
        //   new_in_hole.end());
        // points_before_densification_.insert(h.begin(), h.end());

        holes_v_dens.emplace_back(std::move(hole_dens));
      }
      const Polygon_with_holes pwh_dens(
        std::move(outer_dens),
        holes_v_dens.begin(),
        holes_v_dens.end());
      gd_dens.push_back(std::move(pwh_dens));
    }
    geodivs_dens.emplace_back(std::move(gd_dens));
  }
  geo_divs_ = std::move(geodivs_dens);
  size_t n_points_aft_densification = n_points();
  std::cerr << "Points After Densification: " << n_points_aft_densification
            << std::endl;
  densification_changes.emplace_back(
    n_points_bef_densification,
    n_points_aft_densification);
  is_simple(__func__);
  timer.stop("Densification (using Delanuay Triangles)");
}
