#include <boost/multi_array.hpp>
#include <iostream>
#include <vector>
#include "interpolate_bilinearly.h"
#include "matrix.h";

#include "project.h"

void project(MapState *map_state)
{
  const unsigned int lx = map_state->lx();
  const unsigned int ly = map_state->ly();
  boost::multi_array<XYPoint, 2> &proj = *map_state->proj();

  // Calculate displacement from proj array
  boost::multi_array<double, 2> xdisp(boost::extents[lx][ly]);
  boost::multi_array<double, 2> ydisp(boost::extents[lx][ly]);
  for (unsigned int i = 0; i < lx; i++) {
    for (unsigned int j=0; j<ly; j++) {
      xdisp[i][j] = proj[i][j].x - i - 0.5;
      ydisp[i][j] = proj[i][j].y - j - 0.5;
    }
  }
  std::vector<GeoDiv> new_geo_divs;
  for (auto gd : map_state->geo_divs()) {

    // For each GeoDiv
    GeoDiv new_gd(gd.id());

    for (auto pwh : gd.polygons_with_holes()) {
      // For each polygon with holes

      Polygon old_ext_ring = pwh.outer_boundary();
      Polygon new_ext_ring;

      for (unsigned int i = 0; i < old_ext_ring.size(); i++) {

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
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); hci++) {
        Polygon old_hole = *hci;
        Polygon new_hole;
        for (unsigned int i = 0; i < old_hole.size(); i++) {

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
  map_state->set_geo_divs(new_geo_divs);
  return;
}

// Applying the function find_graticule to any cartogram point would
// return the coordinates (x, y) of the bottom left point of the original
// (untransformed) square graticule the cartogram point was in. This square
// graticule would always have four vertices of coordinates (x, y), (x + 1, y),
// (x + 1, y + 1) and (x, y + 1).

void project_graticule_centroids(MapState *map_state){

  const unsigned int lx = map_state->lx();
  const unsigned int ly = map_state->ly();

  boost::multi_array<XYPoint, 2> &proj = *map_state->proj();
  boost::multi_array<XYPoint, 2> &graticule_centroids = *map_state->graticule_centroids();

  // Resize multi array if running for the first time
  if (map_state->n_finished_integrations() == 0) {
    graticule_centroids.resize(boost::extents[lx][ly]);
  }

  for (unsigned int i = 0; i < lx - 1; i++){
    for (unsigned int j = 0; j < ly - 1; j++){
      graticule_centroids[i][j].x = i;
      graticule_centroids[i][j].y = j;
    }
  }

  // Project graticule centroids
  for (unsigned int i = 0; i < lx - 1; i++){
    for (unsigned int j = 0; j < ly - 1; j++){
      double x_disp =
        proj[i][j].x - i - 0.5 +
        proj[i + 1][j].x - i - 1.5 +
        proj[i][j + 1].x - i - 0.5 +
        proj[i + 1][j + 1].x - i - 1.5;
      double y_disp =
        proj[i][j].y - j - 0.5 +
        proj[i + 1][j].y - j - 0.5 +
        proj[i][j + 1].y - j - 1.5 +
        proj[i + 1][j + 1].y - j - 1.5;
      x_disp /= 4;
      y_disp /= 4;
      graticule_centroids[i][j].x += x_disp;
      graticule_centroids[i][j].y += y_disp;
    }
  }
}

std::vector<XYPoint> find_triangle(const int x,
                                   const int y,
                                   const int lx,
                                   const int ly)
{

  if (x < 0 || x > lx || y < 0 || y > ly) {
    std::cerr << "ERROR: coordinate outside bounding box in "
              << "find_triangle().\n";
    std::cerr << "x=" << x << ", y=" << y << std::endl;
    exit(1);
  }

  std::vector<XYPoint> triangle_coordinates;

  // Get graticule coordinates and centroid.


  int x0 = floor(x + 0.5) - 0.5;
  int y0 = floor(y + 0.5) - 0.5;

  XYPoint centroid;
  centroid.x = x0 + 0.5;
  centroid.y = y0 + 0.5;

  boost::multi_array<double, 2> vx (boost::extents[2][2]);
  boost::multi_array<double, 2> vy (boost::extents[2][2]);

  for (int i = 0; i < 2; i++){
    for (int j = 0; j < 2; j++){
      vx[i][j] = x0 + i;
      vy[i][j] = y0 + j;
    }
  }

  int c0x = 0;
  int c0y = 0;

  if (x > centroid.x){
    c0x = 1;
  }
  if (y > centroid.y){
    c0y = 1;
  }

  double a = (vy[c0x][c0y] - centroid.y) / (vx[c0x][c0y] - centroid.x);
  double b = centroid.y - a * centroid.x;

  int c1x;
  int c1y;

  if (y >= a * x + b){
    if (c0y == 0){
      c1y = 1;
      c1x = c0x;
    } else {
      c1x = 1 - c0x;
      c1y = c0y;
    }
  } else {
    if (c0y == 0){
      c1x = 1 - c0x;
      c1y = c0y;
    } else {
      c1y = 0;
      c1x = c0x;
    }
  }

  triangle_coordinates.push_back(centroid);

  XYPoint c0;
  c0.x = vx[c0x][c0y];
  c0.y = vy[c0x][c0y];
  triangle_coordinates.push_back(c0);

  XYPoint c1;
  c1.x = vx[c1x][c1y];
  c1.y = vy[c1x][c1y];
  triangle_coordinates.push_back(c1);

  return triangle_coordinates;
}

XYPoint affine_trans(std::vector<XYPoint> *tri,
                     std::vector<XYPoint> *org_tri,
                     double x, double y){

  // For each point, we make the following transformation.
  // Suppose we find that, before the cartogram transformation, a point (x, y)
  // is in the triangle (a, b, c). We want to find its position in
  // the projected triangle (p, q, r). We locally approximate the cartogram
  // transformation by an affine transformation T such that T(a) = p,
  // T(b) = q and T(c) = r. We can think of T as a 3x3 matrix
  //  /t11 t12 t13\
  // | t21 t22 t23 |  such that
  //  \ 0   0   1 /
  //  /t11 t12 t13\   /a1 b1 c1\     /p1 q1 r1\
  // | t21 t22 t23 | | a2 b2 c2 | = | p2 q2 r2 | or TA = P. Hence T = PA^{-1}
  //  \ 0   0   1 /   \ 1  1  1/     \ 1  1  1/
  //                              /b2-c2 c1-b1 b1*c2-b2*c1\
  // We have A^{-1} = (1/det(A)) | c2-a2 a1-c1 a2*c1-a1*c2 |. By multiplying
  //                              \a2-b2 b1-a1 a1*b2-a2*b1/
  // PA^{-1} we obtain t11, t12, t13, t21, t22, t23. The postimage of (x, y) i
  // the unprojected map is then "pre" with coordinates
  // post.x = t11*x + t12*y + t13, pre.y = t21*x + t22*y + t23.

  XYPoint pre;
  pre.x = x;
  pre.y = y;

  // Old triangle (a, b, c) as a matrix, explained earlier as Matrix A
  Matrix abc_mA = ((*org_tri)[0], (*org_tri)[1], (*org_tri)[2]);

  // New triangle (p, q, r) as a matrix, explained earlier as Matrix P
  Matrix pqr_mP = ((*tri)[0], (*tri)[1], (*tri)[2]);

  // Calculating transformation matrix
  Matrix mT = pqr_mP.multiply(abc_mA.inverse());

  // Transforming point and pushing back to temporary_ext_boundary
  XYPoint post = mT.transform_XYPoint(pre);

  return post;
}

void project_with_triangulation(MapState *map_state)
{
  const unsigned int lx = map_state->lx();
  const unsigned int ly = map_state->ly();
  boost::multi_array<XYPoint, 2> &proj = *map_state->proj();
  boost::multi_array<XYPoint, 2> &graticule_centroids = *map_state->graticule_centroids();

  std::vector<GeoDiv> new_geo_divs;

  for (auto gd : map_state->geo_divs()) {

    // For each GeoDiv
    GeoDiv new_gd(gd.id());

    for (auto pwh : gd.polygons_with_holes()) {
      // For each polygon with holes

      Polygon old_ext_ring = pwh.outer_boundary();
      Polygon new_ext_ring;

      for (unsigned int i = 0; i < old_ext_ring.size(); i++) {

        // Update exterior ring coordinates

        std::vector<XYPoint> ext_ring_triangle =
          find_triangle(old_ext_ring[i][0], old_ext_ring[i][1],
                        lx, ly);

        std::vector<XYPoint> transformed_ext_ring_triangle;

        XYPoint centroid;
        centroid.x =
          graticule_centroids[int(ext_ring_triangle[0].x) - 1][int(ext_ring_triangle[0].y) - 1].x;
        centroid.y =
          graticule_centroids[int(ext_ring_triangle[0].x) - 1][int(ext_ring_triangle[0].y) - 1].y;

        XYPoint v1;
        v1.x =
          proj[int(ext_ring_triangle[1].x)][int(ext_ring_triangle[1].y)].x;
        v1.y =
          proj[int(ext_ring_triangle[1].x)][int(ext_ring_triangle[1].y)].y;

        XYPoint v2;
        v2.x =
          proj[int(ext_ring_triangle[2].x)][int(ext_ring_triangle[2].y)].x;
        v2.y =
          proj[int(ext_ring_triangle[2].x)][int(ext_ring_triangle[2].y)].y;

        transformed_ext_ring_triangle.push_back(centroid);
        transformed_ext_ring_triangle.push_back(v1);
        transformed_ext_ring_triangle.push_back(v2);

        XYPoint old_ext_ring_intp = affine_trans(&transformed_ext_ring_triangle,
                                             &ext_ring_triangle,
                                             old_ext_ring[i][0], old_ext_ring[i][1]);

        new_ext_ring.push_back(Point(old_ext_ring_intp.x,
                                     old_ext_ring_intp.y));
      }
      std::vector<Polygon> hole_v;
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); hci++) {
        Polygon old_hole = *hci;
        Polygon new_hole;
        for (unsigned int i = 0; i < old_hole.size(); i++) {

          std::vector<XYPoint> hole_triangle =
            find_triangle(old_hole[i][0], old_hole[i][1],
                          lx, ly);

          std::vector<XYPoint> transformed_hole_triangle;

          XYPoint centroid;
          centroid.x =
            graticule_centroids[int(hole_triangle[0].x) - 1][int(hole_triangle[0].y) - 1].x;
          centroid.y =
            graticule_centroids[int(hole_triangle[0].x) - 1][int(hole_triangle[0].y) - 1].y;

          XYPoint v1;
          v1.x =
            proj[int(hole_triangle[1].x)][int(hole_triangle[1].y)].x;
          v1.y =
            proj[int(hole_triangle[1].x)][int(hole_triangle[1].y)].y;

          XYPoint v2;
          v2.x =
            proj[int(hole_triangle[2].x)][int(hole_triangle[2].y)].x;
          v2.y =
            proj[int(hole_triangle[2].x)][int(hole_triangle[2].y)].y;

          transformed_hole_triangle.push_back(centroid);
          transformed_hole_triangle.push_back(v1);
          transformed_hole_triangle.push_back(v2);

          XYPoint old_hole_intp = affine_trans(&transformed_hole_triangle,
                                               &hole_triangle,
                                               old_hole[i][0], old_hole[i][1]);

          new_hole.push_back(Point(old_hole_intp.x,
                                   old_hole_intp.y));
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
  map_state->set_geo_divs(new_geo_divs);
  return;
}

/*
XYPoint affine_trans(double ax, double bx, double cx,
                     double ay, double by, double cy,
                     double px, double qx, double rx,
                     double py, double qy, double ry,
                     double x, double y){
  XYPoint trans_point;
  return trans_point;
}
*/
/*

// After running the project_graticule function, calling graticule_points[x][y].x
// and graticule_points[x][y].y would return the transformed coordinates tx and ty
// of the transformed graticule point that was originally of coordinates (x, y).

void project_graticule(MapState *map_state)
{
  const unsigned int lx = map_state->lx();
  const unsigned int ly = map_state->ly();
  boost::multi_array<XYPoint, 2> &proj = *map_state->proj();

  boost::multi_array<XYPoint, 2> &graticule_points = *map_state->graticule_points();
  boost::multi_array<XYPoint, 2> &graticule_centroids = *map_state->graticule_centroids();

  // Resize multi array if running for the first time
  if (map_state->n_finished_integrations() == 0) {
    graticule_points.resize(boost::extents[lx + 1][ly + 1]);
    graticule_centroids.resize(boost::extents[lx][ly]);
  }

  for (unsigned int i = 0; i < lx + 1; i++){
    for (unsigned int j = 0; j < ly + 1; j++){
      graticule_points[i][j].x = i;
      graticule_points[i][j].y = j;
    }
  }

  for (unsigned int i = 0; i < lx; i++){
    for (unsigned int j = 0; j < ly; j++){
      graticule_centroids[i][j].x = i + 0.5;
      graticule_centroids[i][j].y = j + 0.5;
    }
  }

  // Calculate displacement from proj array
  boost::multi_array<double, 2> xdisp(boost::extents[lx][ly]);
  boost::multi_array<double, 2> ydisp(boost::extents[lx][ly]);
  for (unsigned int i = 0; i < lx; i++) {
    for (unsigned int j=0; j<ly; j++) {
      xdisp[i][j] = proj[i][j].x - i - 0.5;
      ydisp[i][j] = proj[i][j].y - j - 0.5;
    }
  }

  // Project graticule points
  for (unsigned int i = 0; i < lx; i++){
    for (unsigned int j = 0; j < ly; j++){
      double new_x =
        interpolate_bilinearly(graticule_points[i][j].x,
                               graticule_points[i][j].y,
                               &xdisp, 'x', lx, ly);
      double new_y =
        interpolate_bilinearly(graticule_points[i][j].x,
                               graticule_points[i][j].y,
                               &ydisp, 'y', lx, ly);
      graticule_points[i][j].x = new_x;
      graticule_points[i][j].y = new_y;
    }
  }

  // Project graticule centroids
  for (unsigned int i = 0; i < lx; i++){
    for (unsigned int j = 0; j < ly; j++){
      double x_disp =
        graticule_points[i][j].x - i +
        graticule_points[i + 1][j].x - i - 1 +
        graticule_points[i][j + 1].x - i +
        graticule_points[i + 1][j + 1].x - i - 1;
      double y_disp =
        graticule_points[i][j].y - j +
        graticule_points[i + 1][j].y - j +
        graticule_points[i][j + 1].y - j - 1 +
        graticule_points[i + 1][j + 1].y - j - 1;
      x_disp /= 4;
      y_disp /= 4;
      graticule_centroids[i][j].x += x_disp;
      graticule_centroids[i][j].y += y_disp;
    }
  }

  return;
}

// Applying the function find_graticule to any cartogram point would
// return the coordinates (x, y) of the bottom left point of the original
// (untransformed) square graticule the cartogram point was in. This square
// graticule would always have four vertices of coordinates (x, y), (x + 1, y),
// (x + 1, y + 1) and (x, y + 1).

void project_with_triangulation(MapState *map_state)
{
  return;
}
*/
