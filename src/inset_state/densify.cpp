#include "inset_state.hpp"
#include "round_point.hpp"
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <variant>

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
  is_simple(__func__);
  timer.stop("Densification (using Delanuay Triangles)");
}
