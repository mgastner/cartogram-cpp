#include "densification_points.h"
#include "interpolate_bilinearly.h"
#include "matrix.h"
#include "project.h"
#include <boost/multi_array.hpp>
#include <iostream>
#include <vector>

void project(InsetState *inset_state)
{
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  boost::multi_array<XYPoint, 2> &proj = *inset_state->ref_to_proj();
  boost::multi_array<XYPoint, 2> &cum_proj = *inset_state->ref_to_cum_proj();

  // Calculate displacement from proj array
  boost::multi_array<double, 2> xdisp(boost::extents[lx][ly]);
  boost::multi_array<double, 2> ydisp(boost::extents[lx][ly]);
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j=0; j<ly; ++j) {
      xdisp[i][j] = proj[i][j].x - i - 0.5;
      ydisp[i][j] = proj[i][j].y - j - 0.5;
    }
  }

  // Cumulative projection
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {

      // Calculate displacement for cumulative graticule coordinates
      double graticule_intp_x = interpolate_bilinearly(cum_proj[i][j].x,
                                                       cum_proj[i][j].y,
                                                       &xdisp, 'x', lx, ly);
      double graticule_intp_y = interpolate_bilinearly(cum_proj[i][j].x,
                                                       cum_proj[i][j].y,
                                                       &ydisp, 'y', lx, ly);

      // Update cumulative graticule coordinates
      cum_proj[i][j].x += graticule_intp_x;
      cum_proj[i][j].y += graticule_intp_y;
    }
  }
  std::vector<GeoDiv> new_geo_divs;
  for (auto gd : inset_state->geo_divs()) {

    // For each GeoDiv
    GeoDiv new_gd(gd.id());
    for (auto pwh : gd.polygons_with_holes()) {

      // For each polygon with holes
      Polygon old_ext_ring = pwh.outer_boundary();
      Polygon new_ext_ring;
      for (unsigned int i = 0; i < old_ext_ring.size(); ++i) {

        // Update exterior ring coordinates
        double old_ext_intp_x =
          interpolate_bilinearly(old_ext_ring[i][0], old_ext_ring[i][1],
                                 &xdisp, 'x',
                                 lx, ly);
        double old_ext_intp_y =
          interpolate_bilinearly(old_ext_ring[i][0], old_ext_ring[i][1],
                                 &ydisp, 'y',
                                 lx, ly);
        new_ext_ring.push_back(Point(old_ext_intp_x +old_ext_ring[i][0],
                                     old_ext_intp_y + old_ext_ring[i][1]));
      }
      std::vector<Polygon> hole_v;
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
        Polygon old_hole = *hci;
        Polygon new_hole;
        for (unsigned int i = 0; i < old_hole.size(); ++i) {

          // Update hole coordinates
          double old_hole_intp_x =
            interpolate_bilinearly(old_hole[i][0], old_hole[i][1],
                                   &xdisp, 'x',
                                   lx, ly);
          double old_hole_intp_y =
            interpolate_bilinearly(old_hole[i][0], old_hole[i][1],
                                   &ydisp, 'y', lx, ly);
          new_hole.push_back(Point(old_hole_intp_x + old_hole[i][0],
                                   old_hole_intp_y + old_hole[i][1]));
        }
        hole_v.push_back(new_hole);
      }
      const Polygon_with_holes new_pwh(new_ext_ring,
                                       hole_v.begin(),
                                       hole_v.end());
      new_gd.push_back(new_pwh);
    }
    new_geo_divs.push_back(new_gd);
  }
  inset_state->set_geo_divs(new_geo_divs);
  return;
}

// In some chosen_diag() and transformed_triangle(), input x-coordinates can
// only be 0, lx, or 0.5, 1.5, ..., lx-0.5. A similar rule applies to the
// y-coordinates.
void exit_if_point_not_on_grid_or_edge(const XYPoint pt,
                                       InsetState *inset_state)
{
  if ((pt.x != 0.0 &&
       pt.x != inset_state->lx() &&
       pt.x - static_cast<int>(pt.x) != 0.5) ||
      (pt.y != 0.0 &&
       pt.y != inset_state->ly() &&
       pt.y - static_cast<int>(pt.y) != 0.5)) {
    std::cerr << "Error: Invalid input coordinate in triangulation\n"
              << "\tpt = ("
              << pt.x
              << ", "
              << pt.y
              << ")"
              << std::endl;
    exit(1);
  }
  return;
}

XYPoint projected_point_on_grid_or_edge(const XYPoint pt,
                                        InsetState *inset_state)
{
  exit_if_point_not_on_grid_or_edge(pt, inset_state);
  boost::multi_array<XYPoint, 2> &proj = *inset_state->ref_to_proj();
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  XYPoint transf_pt;
  const unsigned int proj_x =
    std::min(static_cast<int>(lx) - 1, static_cast<int>(pt.x));
  const unsigned int proj_y =
    std::min(static_cast<int>(ly) - 1, static_cast<int>(pt.y));
  transf_pt.x = (pt.x == 0.0 || pt.x == inset_state->lx()) ?
                pt.x :
                proj[proj_x][proj_y].x;
  transf_pt.y = (pt.y == 0.0 || pt.y == inset_state->ly()) ?
                pt.y :
                proj[proj_x][proj_y].y;
  return transf_pt;
}

// TODO: chosen_diag() seems to be more naturally thought of as a boolean
// than an integer.

// For a graticule cell with corners stored in the XYPoint array v, determine
// whether the diagonal from v[0] to v[2] is inside the graticule cell. If
// yes, return 0. Otherwise, if the diagonal from v[1] to v[3] is inside the
// graticule cell, return 1. If neither of the two diagonals is inside the
// graticule cell, then the cell's topology is invalid; thus, we exit with an
// error message.
int chosen_diag(const XYPoint v[4],
                unsigned int *n_concave,
                InsetState *inset_state)
{
  // Transform the coordinates in v to the corresponding coordinates on the
  // projected grid. If the x-coordinate is 0 or lx, we keep the input. The
  // input v[i].x can only be 0, lx, or 0.5, 1.5, ..., lx-0.5. A similar rule
  // applies to the y-coordinates. This condition is checked in
  // projected_point_on_grid_or_edge().
  XYPoint tv[4];
  for (unsigned int i = 0; i < 4; ++i) {
    tv[i] = projected_point_on_grid_or_edge(v[i], inset_state);
  }

  // Get the two possible midpoints
  XYPoint midpoint0;
  midpoint0.x = (tv[0].x + tv[2].x) / 2;
  midpoint0.y = (tv[0].y + tv[2].y) / 2;
  XYPoint midpoint1;
  midpoint1.x = (tv[1].x + tv[3].x) / 2;
  midpoint1.y = (tv[1].y + tv[3].y) / 2;

  // Get the transformed graticule cell as a polygon
  Polygon trans_graticule;
  for (unsigned int i = 0; i < 4; i++) {
    trans_graticule.push_back(Point(tv[i].x, tv[i].y));
  }

  // Check if graticule cell is concave
  if (!trans_graticule.is_convex()) {
    ++n_concave;
  }
  if (trans_graticule.bounded_side(Point(midpoint0.x, midpoint0.y)) ==
      CGAL::ON_BOUNDED_SIDE) {
    return 0;
  }
  if (trans_graticule.bounded_side(Point(midpoint1.x, midpoint1.y)) ==
      CGAL::ON_BOUNDED_SIDE) {
    return 1;
  }
  std::cerr << "Invalid graticule cell! At\n";
  std::cerr << "(" << tv[0].x << ", " << tv[0].y << ")\n";
  std::cerr << "(" << tv[1].x << ", " << tv[1].y << ")\n";
  std::cerr << "(" << tv[2].x << ", " << tv[2].y << ")\n";
  std::cerr << "(" << tv[3].x << ", " << tv[3].y << ")\n";
  std::cerr << "i: "
            << static_cast<int>(v[0].x)
            << ", j: "
            << static_cast<int>(v[0].y)
            << std::endl;
  exit(1);
}

void fill_graticule_diagonals(InsetState *inset_state)
{
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  boost::multi_array<int, 2> &graticule_diagonals =
    *inset_state->ref_to_graticule_diagonals();

  // Initialize array if running for the first time
  if (graticule_diagonals.shape()[0] != lx ||
      graticule_diagonals.shape()[1] != ly) {
    graticule_diagonals.resize(boost::extents[lx - 1][ly - 1]);
  }
  unsigned int n_concave = 0;  // Count concave graticule cells
  for (unsigned int i = 0; i < lx - 1; ++i) {
    for (unsigned int j = 0; j < ly - 1; ++j) {
      XYPoint v[4];
      v[0].x = i + 0.5;
      v[0].y = j + 0.5;
      v[1].x = i + 1.5;
      v[1].y = j + 0.5;
      v[2].x = i + 1.5;
      v[2].y = j + 1.5;
      v[3].x = i + 0.5;
      v[3].y = j + 1.5;
      graticule_diagonals[i][j] = chosen_diag(v, &n_concave, inset_state);
    }
  }
  std::cerr << "Number of concave graticule cells: "
            << n_concave
            << std::endl;
  return;
}

// TODO: Using an std::vector seems overkill because we know the size of the
// vector. Should we implement this function with C arrays or std::array
// instead?
std::vector<XYPoint> transformed_triangle(const std::vector<XYPoint> tri,
                                          InsetState *inset_state)
{
  std::vector<XYPoint> transf_tri;
  for (const auto &pt : tri) {
    exit_if_point_not_on_grid_or_edge(pt, inset_state);
    XYPoint transf_pt = projected_point_on_grid_or_edge(pt, inset_state);
    transf_tri.push_back(transf_pt);
  }
  return transf_tri;
}

std::vector<XYPoint> triangle_that_contains_point(const Point pt,
                                                  InsetState *inset_state)
{
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  if (pt.x() < 0 || pt.x() > lx || pt.y() < 0 || pt.y() > ly) {
    std::cerr << "ERROR: coordinate outside bounding box in "
              << "triangle_that_contains_point().\n";
    std::cerr << "pt = " << pt << std::endl;
    exit(1);
  }
  boost::multi_array<int, 2> &graticule_diagonals =
    *inset_state->ref_to_graticule_diagonals();

  // Get original graticule coordinates
  XYPoint v[4];
  v[0].x = std::max(0.0, floor(pt.x() + 0.5) - 0.5);
  v[0].y = std::max(0.0, floor(pt.y() + 0.5) - 0.5);
  v[1].x = std::min(static_cast<double>(lx), floor(pt.x() + 0.5) + 0.5);
  v[1].y = v[0].y;
  v[2].x = v[1].x;
  v[2].y = std::min(static_cast<double>(ly), floor(pt.y() + 0.5) + 0.5);
  v[3].x = v[0].x;
  v[3].y = v[2].y;

  // Assuming that the transformed graticule does not have self-intersections,
  // at least one of the diagonals must be completely inside the graticule.
  // We use that diagonal to split the graticule into two triangles.
  int diag;
  if (v[0].x == 0.0 || v[0].y == 0.0 || v[2].x == lx || v[2].y == ly) {
    unsigned int concave = 0;  // Placeholder value
    diag = chosen_diag(v, &concave, inset_state);
  } else {
    diag =
      graticule_diagonals[static_cast<int>(v[0].x)][static_cast<int>(v[0].y)];
  }
  Polygon tri1;
  Polygon tri2;
  if (diag == 0) {
    tri1.push_back(Point(v[0].x, v[0].y));
    tri1.push_back(Point(v[1].x, v[1].y));
    tri1.push_back(Point(v[2].x, v[2].y));
    tri2.push_back(Point(v[0].x, v[0].y));
    tri2.push_back(Point(v[2].x, v[2].y));
    tri2.push_back(Point(v[3].x, v[3].y));
  } else {
    tri1.push_back(Point(v[0].x, v[0].y));
    tri1.push_back(Point(v[1].x, v[1].y));
    tri1.push_back(Point(v[3].x, v[3].y));
    tri2.push_back(Point(v[1].x, v[1].y));
    tri2.push_back(Point(v[2].x, v[2].y));
    tri2.push_back(Point(v[3].x, v[3].y));
  }

  // Determine which untransformed triangle the given point is in.
  // If the point is in neither, an error is raised.
  std::vector<XYPoint> tri_coords;
  if ((tri1.bounded_side(pt) == CGAL::ON_BOUNDED_SIDE) ||
      (tri1.bounded_side(pt) == CGAL::ON_BOUNDARY)) {
    for (unsigned int i = 0; i < tri1.size(); ++i) {
      XYPoint tri_pt;
      tri_pt.x = tri1[i][0];
      tri_pt.y = tri1[i][1];
      tri_coords.push_back(tri_pt);
    }
  } else if ((tri2.bounded_side(pt) == CGAL::ON_BOUNDED_SIDE) ||
             (tri2.bounded_side(pt) == CGAL::ON_BOUNDARY)) {
    for (unsigned int i = 0; i < tri2.size(); ++i) {
      XYPoint tri_pt;
      tri_pt.x = tri2[i][0];
      tri_pt.y = tri2[i][1];
      tri_coords.push_back(tri_pt);
    }
  } else {
    std::cerr << "Point not in graticule cell!" << std::endl;
    exit(1);
  }
  return tri_coords;
}

XYPoint affine_trans(const std::vector<XYPoint> *tri,
                     const std::vector<XYPoint> *org_tri,
                     const Point pt)
{
  // For each point, we make the following transformation. Suppose we find
  // that, before the cartogram transformation, a point (x, y) is in the
  // triangle (a, b, c). We want to find its position in the projected
  // triangle (p, q, r). We locally approximate the cartogram transformation
  // by an affine transformation T such that T(a) = p, T(b) = q and T(c) = r.
  // We can think of T as a 3x3 matrix
  //    -----------
  //   |t11 t12 t13|
  //   |t21 t22 t23|  such that
  //   | 0   0   1 |
  //    -----------
  //    -----------   ----------     ----------
  //   |t11 t12 t13| | a1 b1 c1 |   | p1 q1 r1 |
  //   |t21 t22 t23| | a2 b2 c2 | = | p2 q2 r2 | or TA = P.
  //   | 0   0   1 | | 1  1  1  |   |  1  1  1 |
  //    -----------   ----------     ----------
  // Hence, T = PA^{-1}.
  //                              -----------------------
  //                             |b2-c2 c1-b1 b1*c2-b2*c1|
  // We have A^{-1} = (1/det(A)) |c2-a2 a1-c1 a2*c1-a1*c2|. By multiplying
  //                             |a2-b2 b1-a1 a1*b2-a2*b1|
  //                              -----------------------
  // PA^{-1} we obtain t11, t12, t13, t21, t22, and t23. If the original
  // coordinates are (x, y) on the unprojected map, then the transformed
  // coordinates are:
  // post.x = t11*x + t12*y + t13, post.y = t21*x + t22*y + t23.
  XYPoint pre;
  pre.x = pt.x();
  pre.y = pt.y();

  // Old triangle (a, b, c) expressed as matrix (A in the comment above)
  Matrix abc_m((*org_tri)[0], (*org_tri)[1], (*org_tri)[2]);

  // New triangle (p, q, r) expressed as matrix (P in the comment above)
  Matrix pqr_m((*tri)[0], (*tri)[1], (*tri)[2]);

  // Transformation matrix T
  Matrix t = pqr_m.multiplied_with(abc_m.inverse());

  // Transformed point
  XYPoint post = t.transformed_XYPoint(pre);
  return post;
}

XYPoint projected_point_with_triangulation(const Point pt,
                                           InsetState *inset_state)
{
  // Get the untransformed triangle the point is in
  std::vector<XYPoint> tri = triangle_that_contains_point(pt, inset_state);

  // Get the coordinates of the transformed triangle
  std::vector<XYPoint> transf_tri = transformed_triangle(tri, inset_state);

  // Compute the projected point from the untransformed and transformed
  // triangles
  return affine_trans(&transf_tri, &tri, pt);
}

void project_with_triangulation(InsetState *inset_state)
{
  // Project GeoDivs
  std::vector<GeoDiv> new_geo_divs;
  for (const auto &gd : inset_state->geo_divs()) {

    // For each GeoDiv
    GeoDiv new_gd(gd.id());
    for (const auto &pwh : gd.polygons_with_holes()) {

      // For each polygon with holes
      const Polygon old_ext_ring = pwh.outer_boundary();
      Polygon new_ext_ring;
      for (unsigned int i = 0; i < old_ext_ring.size(); ++i) {

        // Update exterior ring coordinates
        const XYPoint new_ext_ring_pt =
          projected_point_with_triangulation(old_ext_ring[i], inset_state);
        new_ext_ring.push_back(Point(new_ext_ring_pt.x, new_ext_ring_pt.y));
      }
      std::vector<Polygon> hole_v;
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
        const Polygon old_hole = *hci;
        Polygon new_hole;
        for (unsigned int i = 0; i < old_hole.size(); ++i) {

          // Update hole coordinates
          const XYPoint new_hole_pt =
            projected_point_with_triangulation(old_hole[i], inset_state);
          new_hole.push_back(Point(new_hole_pt.x, new_hole_pt.y));
        }
        hole_v.push_back(new_hole);
      }
      const Polygon_with_holes new_pwh(new_ext_ring,
                                       hole_v.begin(),
                                       hole_v.end());
      new_gd.push_back(new_pwh);
    }
    new_geo_divs.push_back(new_gd);
  }
  inset_state->set_geo_divs(new_geo_divs);

  // Cumulative projection
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  boost::multi_array<XYPoint, 2> &cum_proj = *inset_state->ref_to_cum_proj();
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      const Point old_cum_proj(cum_proj[i][j].x, cum_proj[i][j].y);
      cum_proj[i][j] =
        projected_point_with_triangulation(old_cum_proj, inset_state);
    }
  }
  return;
}
