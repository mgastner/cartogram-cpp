#include "interpolate_bilinearly.h"
#include "matrix.h"
#include "round_point.h"
#include <boost/multi_array.hpp>
#include <iostream>
#include <vector>

Point interpolate_point_bilinearly(Point p1,
                                   const boost::multi_array<double, 2> *xdisp,
                                   const boost::multi_array<double, 2> *ydisp,
                                   const unsigned int lx,
                                   const unsigned int ly)
{
  const double intp_x =
    interpolate_bilinearly(p1.x(), p1.y(), xdisp, 'x', lx, ly);
  const double intp_y =
    interpolate_bilinearly(p1.x(), p1.y(), ydisp, 'y', lx, ly);
  return Point(p1.x() + intp_x, p1.y() + intp_y);
}

void InsetState::project()
{
  // Calculate displacement from proj array
  boost::multi_array<double, 2> xdisp(boost::extents[lx_][ly_]);
  boost::multi_array<double, 2> ydisp(boost::extents[lx_][ly_]);
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j=0; j<ly_; ++j) {
      xdisp[i][j] = proj_[i][j].x - i - 0.5;
      ydisp[i][j] = proj_[i][j].y - j - 0.5;
    }
  }

  // Cumulative projection
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {

      // TODO: Should the interpolation be made on the basis of triangulation?
      // Calculate displacement for cumulative graticule coordinates
      const double graticule_intp_x = interpolate_bilinearly(
        cum_proj_[i][j].x,
        cum_proj_[i][j].y,
        &xdisp, 'x', lx_, ly_);
      const double graticule_intp_y = interpolate_bilinearly(
        cum_proj_[i][j].x,
        cum_proj_[i][j].y,
        &ydisp, 'y', lx_, ly_);

      // Update cumulative graticule coordinates
      cum_proj_[i][j].x += graticule_intp_x;
      cum_proj_[i][j].y += graticule_intp_y;
    }
  }

  // Specialise interpolate_point_bilinearly
  std::function<Point(Point)> lambda =
    [&xdisp, &ydisp, lx = lx_, ly = ly_](Point p1) {
      return interpolate_point_bilinearly(p1, &xdisp, &ydisp, lx, ly);
    };

  // Apply "lambda" to all points
  transform_points(lambda);
  return;
}

// In chosen_diag() and transformed_triangle(), the input x-coordinates can
// only be 0, lx, or 0.5, 1.5, ..., lx-0.5. A similar rule applies to the
// y-coordinates.
void exit_if_not_on_grid_or_edge(const Point pt,
                                 const unsigned int lx,
                                 const unsigned int ly)
{
  if ((pt.x() != 0.0 && pt.x() != lx && pt.x() - int(pt.x()) != 0.5) ||
      (pt.y() != 0.0 && pt.y() != ly && pt.y() - int(pt.y()) != 0.5)) {
    std::cerr << "Error: Invalid input coordinate in triangulation\n"
              << "\tpt = ("
              << pt.x()
              << ", "
              << pt.y()
              << ")"
              << std::endl;
    exit(1);
  }
  return;
}

Point InsetState::projected_point(const Point p1)
{
  exit_if_not_on_grid_or_edge(p1, lx_, ly_);
  const unsigned int proj_x = std::min(
    static_cast<unsigned int>(lx_) - 1,
    static_cast<unsigned int>(p1.x()));
  const unsigned int proj_y = std::min(
    static_cast<unsigned int>(ly_) - 1,
    static_cast<unsigned int>(p1.y()));
  return Point((p1.x() == 0.0 || p1.x() == lx_) ?
               p1.x() :
               proj_[proj_x][proj_y].x,
               (p1.y() == 0.0 || p1.y() == ly_) ?
               p1.y() :
               proj_[proj_x][proj_y].y);
}

// TODO: chosen_diag() seems to be more naturally thought of as a boolean
// than an integer.

// For a graticule cell with corners stored in the XYPoint array v, determine
// whether the diagonal from v[0] to v[2] is inside the graticule cell. If
// yes, return 0. Otherwise, if the diagonal from v[1] to v[3] is inside the
// graticule cell, return 1. If neither of the two diagonals is inside the
// graticule cell, then the cell's topology is invalid; thus, we exit with an
// error message.
int InsetState::chosen_diag(const Point v[4], unsigned int *num_concave)
{
  // The input v[i].x can only be 0, lx, or 0.5, 1.5, ..., lx-0.5. A similar
  // rule applies to the y-coordinates.
  for (unsigned int i = 0; i < 4; ++i) {
    exit_if_not_on_grid_or_edge(v[i], lx_, ly_);
  }

  // Transform the coordinates in v to the corresponding coordinates on the
  // projected grid. If the x-coordinate is 0 or lx, we keep the input. The
  // input v[i].x can only be 0, lx, or 0.5, 1.5, ..., lx-0.5. A similar rule
  // applies to the y-coordinates. This condition is checked in
  // projected_point().
  Point tv[4];
  for (unsigned int i = 0; i < 4; ++i) {
    tv[i] = projected_point(v[i]);
  }

  // Get the two possible midpoints
  const Point midpoint0(
    (tv[0].x() + tv[2].x()) / 2,
    (tv[0].y() + tv[2].y()) / 2);
  const Point midpoint1(
    (tv[1].x() + tv[3].x()) / 2,
    (tv[1].y() + tv[3].y()) / 2);

  // Get the transformed graticule cell as a polygon
  Polygon trans_graticule;
  for (unsigned int i = 0; i < 4; ++i) {
    trans_graticule.push_back(tv[i]);
  }

  // Check if graticule cell is concave
  if (!trans_graticule.is_convex()) {
    *num_concave += 1;
  }
  if (trans_graticule.bounded_side(midpoint0) == CGAL::ON_BOUNDED_SIDE) {
    return 0;
  }
  if (trans_graticule.bounded_side(midpoint1) == CGAL::ON_BOUNDED_SIDE) {
    return 1;
  }
  std::cerr << "Invalid graticule cell! At\n";
  std::cerr << "(" << tv[0].x() << ", " << tv[0].y() << ")\n";
  std::cerr << "(" << tv[1].x() << ", " << tv[1].y() << ")\n";
  std::cerr << "(" << tv[2].x() << ", " << tv[2].y() << ")\n";
  std::cerr << "(" << tv[3].x() << ", " << tv[3].y() << ")\n";
  std::cerr << "Original: \n";
  std::cerr << "(" << v[0].x() << ", " << v[0].y() << ")\n";
  std::cerr << "(" << v[1].x() << ", " << v[1].y() << ")\n";
  std::cerr << "(" << v[2].x() << ", " << v[2].y() << ")\n";
  std::cerr << "(" << v[3].x() << ", " << v[3].y() << ")\n";
  std::cerr << "i: "
            << static_cast<unsigned int>(v[0].x())
            << ", j: "
            << static_cast<unsigned int>(v[0].y())
            << std::endl;
  exit(1);
}

void InsetState::fill_graticule_diagonals()
{
  // Initialize array if running for the first time
  if (graticule_diagonals_.shape()[0] != lx_ ||
      graticule_diagonals_.shape()[1] != ly_) {
    graticule_diagonals_.resize(boost::extents[lx_ - 1][ly_ - 1]);
  }
  unsigned int n_concave = 0;  // Count concave graticule cells
  for (unsigned int i = 0; i < lx_ - 1; ++i) {
    for (unsigned int j = 0; j < ly_ - 1; ++j) {
      Point v[4];
      v[0] = Point(double(i) + 0.5, double(j) + 0.5);
      v[1] = Point(double(i) + 1.5, double(j) + 0.5);
      v[2] = Point(double(i) + 1.5, double(j) + 1.5);
      v[3] = Point(double(i) + 0.5, double(j) + 1.5);
      graticule_diagonals_[i][j] = chosen_diag(v, &n_concave);
    }
  }
  std::cerr << "Number of concave graticule cells: "
            << n_concave
            << std::endl;
  return;
}

std::array<Point, 3> InsetState::transformed_triangle(const std::array<Point, 3>
                                                      tri)
{
  std::array<Point, 3> transf_tri;
  for (unsigned int i = 0; i < 3; ++i) {
    exit_if_not_on_grid_or_edge(tri[i], lx_, ly_);
    const auto transf_pt = projected_point(tri[i]);
    transf_tri[i] = transf_pt;
  }
  return transf_tri;
}

// Determine if a point `pt` is on the boundary of a triangle by using cross
// products to find areas spanned by pt and each triangle edge. Idea from:
// https://stackoverflow.com/questions/7050186/find-if-point-lies-on-line-segment
// This function is needed because, sometimes,
// `triangle.bounded_side(Point(x, y)) == CGAL::ON_BOUNDARY` does not return
// `true` even if the point is on the boundary.
bool is_on_triangle_boundary(const Point pt, const Polygon triangle)
{
  for (unsigned int i = 0; i < triangle.size(); ++i) {
    const auto t1 = triangle[i];
    const auto t2 = triangle[(i == triangle.size() - 1) ? 0 : i + 1];
    const double area = (t1.x() - pt.x()) * (t2.y() - pt.y()) -
                        (t2.x() - pt.x()) * (t1.y() - pt.y());
    if (almost_equal(area, 0.0)) {
      return true;
    }
  }
  return false;
}

// Get the untransformed coordinates of the triangle in which the point `pt`
// is located. After transformation, this triangle must be entirely inside
// the transformed graticule cell.
std::array<Point, 3> InsetState::untransformed_triangle(const Point pt)
{
  if (pt.x() < 0 || pt.x() > lx_ || pt.y() < 0 || pt.y() > ly_) {
    CGAL::set_pretty_mode(std::cerr);
    std::cerr << "ERROR: coordinate outside bounding box in "
              << __func__
              << "().\npt = "
              << pt
              << std::endl;
    exit(1);
  }

  // Get original graticule coordinates
  Point v[4];
  v[0] = Point(std::max(0.0, floor(pt.x() + 0.5) - 0.5),
               std::max(0.0, floor(pt.y() + 0.5) - 0.5));
  v[1] = Point(std::min(static_cast<double>(lx_), floor(pt.x() + 0.5) + 0.5),
               v[0].y());
  v[2] = Point(v[1].x(),
               std::min(static_cast<double>(ly_), floor(pt.y() + 0.5) + 0.5));
  v[3] = Point(v[0].x(),
               v[2].y());

  // TODO: diag SEEMS TO BE MORE NATURALLY THOUGHT OF AS bool INSTEAD OF int.
  // Assuming that the transformed graticule does not have self-intersections,
  // at least one of the diagonals must be completely inside the graticule.
  // We use that diagonal to split the graticule into two triangles.
  int diag;
  if (v[0].x() == 0.0 || v[0].y() == 0.0 || v[2].x() == lx_ || v[2].y() == ly_) {

    // Case where the graticule is on the edge of the grid.
    // We calculate the chosen diagonal, as graticule_diagonals does not store
    // the diagonals for edge grid cells.
    unsigned int n_concave = 0;
    diag = chosen_diag(v, &n_concave);
  } else {

    // Case where the graticule is not on the edge of the grid. We can find the
    // already computed chosen diagonal in graticule_diagonals_.
    const unsigned int x = static_cast<unsigned int>(v[0].x());
    const unsigned int y = static_cast<unsigned int>(v[0].y());
    diag = graticule_diagonals_[x][y];
  }

  // Get the two possible triangles
  Polygon triangle1;
  Polygon triangle2;
  if (diag == 0) {
    triangle1.push_back(v[0]);
    triangle1.push_back(v[1]);
    triangle1.push_back(v[2]);
    triangle2.push_back(v[0]);
    triangle2.push_back(v[2]);
    triangle2.push_back(v[3]);
  } else {
    triangle1.push_back(v[0]);
    triangle1.push_back(v[1]);
    triangle1.push_back(v[3]);
    triangle2.push_back(v[1]);
    triangle2.push_back(v[2]);
    triangle2.push_back(v[3]);
  }

  // Determine which untransformed triangle the given point is in. If the
  // point is in neither, an error is raised.
  std::array<Point, 3> triangle_coordinates;
  if ((triangle1.bounded_side(pt) == CGAL::ON_BOUNDED_SIDE) ||
      is_on_triangle_boundary(pt, triangle1)) {
    for (unsigned int i = 0; i < triangle1.size(); ++i) {
      triangle_coordinates[i] = triangle1[i];
    }
  } else if ((triangle2.bounded_side(pt) == CGAL::ON_BOUNDED_SIDE) ||
             is_on_triangle_boundary(pt, triangle2)) {
    for (unsigned int i = 0; i < triangle2.size(); ++i) {
      triangle_coordinates[i] = triangle2[i];
    }
  } else {
    std::cerr << "Point not in graticule cell!\n";
    std::cerr << "Point coordinates:\n";
    std::cerr << "(" << pt.x() << ", " << pt.y() << ")\n";
    std::cerr << "Original graticule cell:\n";
    std::cerr << "(" << v[0].x() << ", " << v[0].y() << ")\n";
    std::cerr << "(" << v[1].x() << ", " << v[1].y() << ")\n";
    std::cerr << "(" << v[2].x() << ", " << v[2].y() << ")\n";
    std::cerr << "(" << v[3].x() << ", " << v[3].y() << ")\n";
    std::cerr << "Chosen diagonal: " << diag << "\n";
    exit(1);
  }
  return triangle_coordinates;
}

Point affine_trans(const std::array<Point, 3> tri,
                   const std::array<Point, 3> org_tri,
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
  const Point pre(pt.x(),pt.y());

  // Old triangle (a, b, c) expressed as matrix A
  const Matrix abc_mA(org_tri[0], org_tri[1], org_tri[2]);

  // New triangle (p, q, r) expressed as matrix P
  const Matrix pqr_mP(tri[0], tri[1], tri[2]);

  // Transformation matrix T
  const auto mT = pqr_mP.multiplied_with(abc_mA.inverse());

  // Transformed point
  return mT.transformed_point(pre);
}

Point InsetState::projected_point_with_triangulation(const Point pt)
{
  // Get the untransformed triangle the point is in
  const auto old_triangle = untransformed_triangle(pt);

  // Get the coordinates of the transformed triangle
  const auto new_triangle = transformed_triangle(old_triangle);

  // Get the transformed point and return it
  const auto transformed_pt = affine_trans(new_triangle, old_triangle, pt);
  return rounded_point(transformed_pt, lx_, ly_);
}

void InsetState::project_with_triangulation()
{

  // Iterate over GeoDivs
  for (auto &gd : geo_divs_) {

    // Iterate over Polygon_with_holes
    for (auto &pwh : *gd.ref_to_polygons_with_holes()) {

      // Get outer boundary
      auto &outer_boundary = *(&pwh.outer_boundary());

      // Iterate over outer boundary's coordinates
      for (auto &coords_outer : outer_boundary) {

        // Assign outer boundary's coordinates to transformed coordinates
        coords_outer = projected_point_with_triangulation(coords_outer);
      }

      // Iterate over holes
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {

        // Iterate over hole's coordinates
        for (auto &coords_hole : *h) {

          // Assign hole's coordinates to transformed coordinates
          coords_hole = projected_point_with_triangulation(coords_hole);
        }
      }
    }
  }

  // Cumulative projection
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      const Point old_cum_proj(cum_proj_[i][j].x, cum_proj_[i][j].y);
      const auto new_cum_proj_pt =
        projected_point_with_triangulation(old_cum_proj);
      cum_proj_[i][j].x = new_cum_proj_pt.x();
      cum_proj_[i][j].y = new_cum_proj_pt.y();
    }
  }
  return;
}
