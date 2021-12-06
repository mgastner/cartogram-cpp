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
  boost::multi_array<XYPoint, 2> &proj = *inset_state->proj();

  // Calculate displacement from proj array
  boost::multi_array<double, 2> xdisp(boost::extents[lx][ly]);
  boost::multi_array<double, 2> ydisp(boost::extents[lx][ly]);
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j=0; j<ly; ++j) {
      xdisp[i][j] = proj[i][j].x - i - 0.5;
      ydisp[i][j] = proj[i][j].y - j - 0.5;
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
void exit_if_point_not_on_grid_or_edge(Point pt, InsetState *inset_state)
{
  if ((pt.x() != 0.0 && pt.x() != inset_state->lx() && pt.x() - int(pt.x()) != 0.5) ||
      (pt.y() != 0.0 && pt.y() != inset_state->ly() && pt.y() - int(pt.y()) != 0.5)) {
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

// For a graticule cell with corners stored in the XYPoint array v, determine
// whether the diagonal from v[0] to v[2] is inside the graticule cell. If
// yes, return 0. Otherwise, if the diagonal from v[1] to v[3] is inside the
// graticule cell, return 1. If neither of the two diagonals is inside the
// graticule cell, then the cell's topology is invalid; thus, we exit with an
// error message.
int chosen_diag(Point v[4], int *num_concave, InsetState *inset_state)
{
  boost::multi_array<XYPoint, 2> &proj = *inset_state->proj();
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();

  // The input v[i].x can only be 0, lx, or 0.5, 1.5, ..., lx-0.5.
  // A similar rule applies to the y-coordinates.
  for (unsigned int i = 0; i < 4; i++) {
    exit_if_point_not_on_grid_or_edge(v[i], inset_state);
  }

  // Transform the coordinates in v to the corresponding coordinates on the
  // projected grid. If the x-coordinate is 0 or lx, we keep the input.
  XYPoint tv[4];
  for (unsigned int i = 0; i < 4; i++) {
    int proj_x = std::min(int(lx) - 1, int(v[i].x()));
    int proj_y = std::min(int(ly) - 1, int(v[i].y()));
    tv[i].x = (v[i].x() != 0 && v[i].x() != lx) ?
              proj[proj_x][proj_y].x :
              v[i].x();
    tv[i].y = (v[i].y() != 0 && v[i].y() != ly) ?
              proj[proj_x][proj_y].y :
              v[i].y();
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
    *num_concave += 1;
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
  std::cerr << "Original: \n";
  std::cerr << "(" << v[0].x() << ", " << v[0].y() << ")\n";
  std::cerr << "(" << v[1].x() << ", " << v[1].y() << ")\n";
  std::cerr << "(" << v[2].x() << ", " << v[2].y() << ")\n";
  std::cerr << "(" << v[3].x() << ", " << v[3].y() << ")\n";
  std::cerr << "Transformed: \n";
  std::cerr << "(" << proj[v[0].x()][v[0].y()].x << ", " << proj[v[0].x()][v[0].y()].y << ")\n";
  std::cerr << "(" << proj[v[1].x()][v[1].y()].x << ", " << proj[v[1].x()][v[1].y()].y << ")\n";
  std::cerr << "(" << proj[v[2].x()][v[2].y()].x << ", " << proj[v[2].x()][v[2].y()].y << ")\n";
  std::cerr << "(" << proj[v[3].x()][v[3].y()].x << ", " << proj[v[3].x()][v[3].y()].y << ")\n";
  std::cerr << "i: " << int(v[0].x()) << "j: " << int(v[0].y()) << "\n";
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
    graticule_diagonals.resize(boost::extents[lx][ly]);
  }
  int num_concave = 0;  // Count concave graticule cells
  for (unsigned int i = 0; i < lx - 1; ++i) {
    for (unsigned int j = 0; j < ly - 1; ++j) {
      Point v[4];
      v[0] = Point(double(i) + 0.5, double(j) + 0.5);
      v[1] = Point(double(i) + 1.5, double(j) + 0.5);
      v[2] = Point(double(i) + 1.5, double(j) + 1.5);
      v[3] = Point(double(i) + 0.5, double(j) + 1.5);
      graticule_diagonals[i][j] = chosen_diag(v, &num_concave, inset_state);
    }
  }
  std::cerr << "Number of concave graticule cells: " << num_concave << "\n";
  return;
}

// TODO: Using an std::vector seems overkill because we know the size of the
// vector. Should we implement this function with C arrays instead?
std::array<Point, 3> transformed_triangle(std::array<Point, 3> triangle,
                                          InsetState *inset_state)
{
  boost::multi_array<XYPoint, 2> &proj = *inset_state->proj();
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  
  std::array<Point, 3> transformed_triangle;
  for(unsigned int i = 0; i < 3; ++i) {
    exit_if_point_not_on_grid_or_edge(triangle[i], inset_state);
    int proj_x = std::min(int(lx) - 1, int(triangle[i].x()));
    int proj_y = std::min(int(ly) - 1, int(triangle[i].y()));
    Point transformed_point(
      (triangle[i].x() == 0.0 || triangle[i].x() == inset_state->lx()) ?
      triangle[i].x() : proj[proj_x][proj_y].x,
      (triangle[i].y() == 0.0 || triangle[i].y() == inset_state->ly()) ?
      triangle[i].y() : proj[proj_x][proj_y].y
    );
    transformed_triangle[i] = transformed_point;
  }
  return transformed_triangle;
}

bool is_point_on_triangle_boundary(Polygon triangle,
                                   const double x,
                                   const double y)
{
  // From: https://stackoverflow.com/questions/7050186/find-if-point-lies-on-line-segment
  for (unsigned int i = 0; i < triangle.size(); i++){
    double tx1 = triangle[i].x();
    double ty1 = triangle[i].y();
    double tx2 = triangle[(i == triangle.size() - 1) ? 0 : i + 1].x();
    double ty2 = triangle[(i == triangle.size() - 1) ? 0 : i + 1].y();
    double seg = std::sqrt((tx1 - tx2) * (tx1 - tx2) + (ty1 - ty2) * (ty1 - ty2));
    double seg1 = std::sqrt((tx1 - x) * (tx1 - x) + (ty1 - y) * (ty1 - y));
    double seg2 = std::sqrt((tx2 - x) * (tx2 - x) + (ty2 - y) * (ty2 - y));
    if (almost_equal(seg, seg1 + seg2)) return true;
  }
  return false;
}

std::array<Point, 3> untransformed_triangle(const double x,
                                   const double y,
                                   InsetState *inset_state)
{
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  if (x < 0 || x > lx || y < 0 || y > ly) {
    std::cerr << "ERROR: coordinate outside bounding box in "
              << "find_triangle().\n";
    std::cerr << "x=" << x << ", y=" << y << std::endl;
    std::cerr << almost_equal(x, 0.0) << "\n";
    exit(1);
  }
  boost::multi_array<int, 2> &graticule_diagonals =
    *inset_state->ref_to_graticule_diagonals();

  // Get original graticule coordinates
  Point v[4];
  v[0] = Point(std::max(0.0, floor(x + 0.5) - 0.5),
               std::max(0.0, floor(y + 0.5) - 0.5));
  v[1] = Point(std::min(double(lx), floor(x + 0.5) + 0.5),
               v[0].y());
  v[2] = Point(std::min(double(lx), floor(x + 0.5) + 0.5),
               std::min(double(ly), floor(y + 0.5) + 0.5));
  v[3] = Point(v[0].x(),
               std::min(double(ly), floor(y + 0.5) + 0.5));

  // Assuming that the transformed graticule does not have self-intersections,
  // at least one of the diagonals must be completely inside the graticule.
  // We use that diagonal to split the graticule into two triangles.
  int diag;
  if (v[0].x() == 0.0 || v[0].y() == 0.0 || v[2].x() == lx || v[2].y() == ly) {
    int concave = 0;
    diag = chosen_diag(v, &concave, inset_state);
  } else {
    diag = graticule_diagonals[int(v[0].x())][int(v[0].y())];
  }
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

  // Determine which untransformed triangle the given point is in.
  // If the point is in neither, an error is raised.
  std::array<Point, 3> triangle_coordinates;
  if ((triangle1.bounded_side(Point(x, y)) == CGAL::ON_BOUNDED_SIDE) ||
      (is_point_on_triangle_boundary(triangle1, x, y))) {
  // if ((triangle1.bounded_side(Point(x, y)) == CGAL::ON_BOUNDED_SIDE) ||
  //     (triangle1.bounded_side(Point(x, y)) == CGAL::ON_BOUNDARY)) {
    for (unsigned int i = 0; i < triangle1.size(); ++i) {
      triangle_coordinates[i] = triangle1[i];
    }
  } else if ((triangle2.bounded_side(Point(x, y)) == CGAL::ON_BOUNDED_SIDE) ||
             (is_point_on_triangle_boundary(triangle2, x, y))) {
    // } else if ((triangle2.bounded_side(Point(x, y)) == CGAL::ON_BOUNDED_SIDE) ||
    //          (triangle2.bounded_side(Point(x, y)) == CGAL::ON_BOUNDARY)) {
    for (unsigned int i = 0; i < triangle2.size(); ++i) {
      triangle_coordinates[i] = triangle2[i];
    }
  } else {
    std::cerr << "Point not in graticule cell!\n";
    std::setprecision(20);
    Point p = Point(x, y);
    std::cerr << "Point coordinates:\n";
    std::cerr << "(" << p.x() << ", " << p.y() << ")\n";
    std::cerr << "Original graticule cell:\n";
    std::cerr << "(" << v[0].x() << ", " << v[0].y() << ")\n";
    std::cerr << "(" << v[1].x() << ", " << v[1].y() << ")\n";
    std::cerr << "(" << v[2].x() << ", " << v[2].y() << ")\n";
    std::cerr << "(" << v[3].x() << ", " << v[3].y() << ")\n";
    std::cerr << "Chosen diagonal: " << diag << "\n";
    std::cerr << "Triangle 1:\n";
    for (unsigned int i = 0; i < triangle1.size(); i++){
      std::cerr << "(" << triangle1[i].x() << ", " << triangle1[i].y() << ")\n";
    }
    std::cerr << (triangle1.bounded_side(Point(x, y)) == CGAL::ON_BOUNDARY) << "\n";
    std::cerr << (triangle1.bounded_side(Point(x, y)) == CGAL::ON_BOUNDED_SIDE) << "\n";
    std::cerr << (triangle2.bounded_side(Point(x, y)) == CGAL::ON_BOUNDARY) << "\n";
    std::cerr << (triangle2.bounded_side(Point(x, y)) == CGAL::ON_BOUNDED_SIDE) << "\n";
    std::cerr << (p.y() == triangle1[0].y()) << "\n";
    std::cerr << almost_equal(p.y(), triangle1[0].y()) << "\n";
    std::cerr << is_point_on_triangle_boundary(triangle1, x, y) << "\n";
    exit(1);
  }
  return triangle_coordinates;
}

Point affine_trans(std::array<Point, 3> tri,
                   std::array<Point, 3> org_tri,
                   Point pt)
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
  // PA^{-1} we obtain t11, t12, t13, t21, t22, t23. If the original
  // coordinates are (x, y) on the unprojected map, then the transformed
  // coordinates are:
  // post.x = t11*x + t12*y + t13, post.y = t21*x + t22*y + t23.
  Point pre(pt.x(),pt.y());

  // Old triangle (a, b, c) expressed as matrix A
  Matrix abc_mA(org_tri[0], org_tri[1], org_tri[2]);

  // New triangle (p, q, r) expressed as matrix P
  Matrix pqr_mP(tri[0], tri[1], tri[2]);

  // Transformation matrix T
  Matrix mT = pqr_mP.multiplied_with(abc_mA.inverse());

  // Transformed point
  Point post = mT.transformed_point(pre);
  return post;
}

Point transformed_point(Point old_point, InsetState *inset_state)
{
  // Get the untransformed triangle the point is in.
  std::array<Point, 3> old_triangle =
    untransformed_triangle(old_point.x(), old_point.y(), inset_state);
  
  // Get the coordinates of the transformed triangle.
  std::array<Point, 3> new_triangle =
    transformed_triangle(old_triangle, inset_state);

  // Get the transformed point and return it.
  Point transformed_point = affine_trans(new_triangle, old_triangle,
                                         old_point);                                     
  return transformed_point;
}

void project_with_triangulation(InsetState *inset_state)
{
  std::vector<GeoDiv> new_geo_divs;
  for (auto gd : inset_state->geo_divs()) {

    // For each GeoDiv
    GeoDiv new_gd(gd.id());
    for (auto pwh : gd.polygons_with_holes()) {

      // For each polygon with holes
      Polygon old_ext_ring = pwh.outer_boundary();
      Polygon new_ext_ring;
      for (unsigned int i = 0; i < old_ext_ring.size(); ++i) {
        new_ext_ring.push_back(transformed_point(old_ext_ring[i],
                                                 inset_state));
      }
      std::vector<Polygon> hole_v;
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
        Polygon old_hole = *hci;
        Polygon new_hole;
        for (unsigned int i = 0; i < old_hole.size(); ++i) {
          new_hole.push_back(transformed_point(old_hole[i], inset_state));
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
