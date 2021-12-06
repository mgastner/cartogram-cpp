#include "densification_points.h"

// Use machine epsilon (defined in constants.h) to get almost equal doubles.
// From https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
bool almost_equal(double a, double b) {
  return fabs(a - b) <= dbl_epsilon * fabs(a + b) * 2;
}

// Determine whether points are indistinguishable
bool points_almost_equal(Point a, Point b) {
  return (almost_equal(a[0], b[0]) && almost_equal(a[1], b[1]));
}
bool xy_points_almost_equal(XYPoint a, XYPoint b) {
  return (almost_equal(a.x, b.x) && almost_equal(a.y, b.y));
}

// This function takes two lines as input:
// - line `a`, defined by points a1 and a2.
// - line `b`, defined by points b1 and b2.
// The function returns the intersection between them. If the two lines are
// parallel or are the same, the function returns the point (-1, -1), which
// is always outside of any graticule grid cell.
XYPoint calc_intersection(XYPoint a1, XYPoint a2, XYPoint b1, XYPoint b2) {

  // Check if any segment is undefined (i.e., defined by identical points)
  // if (a1 == a2 || b1 == b2) {
  //   std::cerr << "ERROR: End points of line segment are identical"
  //             << std::endl;
  //   _Exit(EXIT_FAILURE);
  // }

  // Get line equations
  double a = (a1.y - a2.y) / (a1.x - a2.x);
  double a_intercept = a1.y - (a1.x * a);
  double b = (b1.y - b2.y) / (b1.x - b2.x);
  double b_intercept = b1.y - (b1.x * b);
  XYPoint intersection;
  if (isfinite(a) && isfinite(b) && a != b) {

    // Neither the line (a1, a2) nor the line (b1, b2) is vertical
    intersection.x = (b_intercept - a_intercept) / (a - b);
    intersection.y = a * intersection.x + a_intercept;
  } else if (isfinite(a) && isinf(b)) {

    // Only line (b1, b2) is vertical
    intersection.x = b1.x;
    intersection.y = a * b1.x + a_intercept;

  } else if (isfinite(b) && isinf(a)) {

    // Only line (a1, a2) is vertical
    intersection.x = a1.x;
    intersection.y = b * a1.x + b_intercept;

  } else {

    // Set negative intersection coordinates if there is no solution or
    // infinitely many solutions
    intersection.x = -1;
    intersection.y = -1;
  }
  return intersection;
}

// TODO: If a_ or b_ are themselves intersection points (e.g. if pt1.x is an
// integer plus 0.5), it appears to be included in the returned intersections.
// Would this property cause the point to be included twice in the line
// segment (once when the end point is the argument pt1 and a second time when
// the same end point is the argument pt2)?

// This function takes two points (called pt1 and pt2) and returns all
// horizontal and vertical intersections of the line segment between pt1 and
// pt2 with a graticule whose graticule lines are placed one unit apart. The
// function also returns all intersections with the diagonals of these
// graticule cells. The function assumes that graticule cells start at
// (0.5, 0.5).
std::vector<Point> densification_points(Point pt1, Point pt2,
                                        const unsigned int lx,
                                        const unsigned int ly)
{
  // Vector for storing intersections before removing duplicates
  std::vector<XYPoint> temp_intersections;

  // Assign `a` as the leftmost point of p1 and pt2. If both points have the
  // same x-coordinate, assign `a` as the lower point. The other point is
  // assigned as b. The segments (a, b) and (b, a) describe the same segment.
  // However, if we flip the order of a and b, the resulting intersections are
  // not necessarily the same because of floating point errors.
  XYPoint a;
  XYPoint b;
  if ((pt1[0] > pt2[0]) || ((pt1[0] == pt2[0]) && (pt1[1] > pt2[1]))) {
    a.x = pt2[0]; a.y = pt2[1];
    b.x = pt1[0]; b.y = pt1[1];
  } else{
    a.x = pt1[0]; a.y = pt1[1];
    b.x = pt2[0]; b.y = pt2[1];
  }
  temp_intersections.push_back(a);
  temp_intersections.push_back(b);

  // Get bottom-left point of graticule cell containing `a`
  XYPoint av0;
  // av0.x = floor(a.x + 0.5) - 0.5;
  // av0.y = floor(a.y + 0.5) - 0.5;
  av0.x = std::max(0.0, floor(a.x + 0.5) - 0.5);
  av0.y = std::max(0.0, floor(a.y + 0.5) - 0.5);

  // Get bottom-left point of graticule cell containing `b`
  XYPoint bv0;
  // bv0.x = floor(b.x + 0.5) - 0.5;
  // bv0.y = floor(b.y + 0.5) - 0.5;
  bv0.x = std::max(0.0, floor(b.x + 0.5) - 0.5);
  bv0.y = std::max(0.0, floor(b.y + 0.5) - 0.5);

  // Get bottom-left (start_v) and top-right (end_v) graticule cells of the
  // graticule cell rectangle (the smallest rectangular section of the
  // graticule grid cell containing both points)
  XYPoint start_v;
  XYPoint end_v;
  start_v.x = av0.x;
  end_v.x = bv0.x;
  if (a.y <= b.y) {
    start_v.y = av0.y;
    end_v.y = bv0.y;
  } else {
    start_v.y = bv0.y;
    end_v.y = av0.y;
  }

  int dist_x = std::ceil(end_v.x - start_v.x);
  int dist_y = std::ceil(end_v.y - start_v.y);
  double cur_x = start_v.x;
  double cur_y = start_v.y;

  // Loop through each row, from bottom to top
  // for (unsigned int int_y = (unsigned int)(start_v.y);
  //      int_y <= (unsigned int)(end_v.y);
  //      ++int_y) {
  for (int i = 0; i <= dist_y; ++i){

    // Loop through each column, from left to right
    // for (unsigned int int_x = (unsigned int)(start_v.x);
    //      int_x <= (unsigned int)(end_v.x);
    //      ++int_x) {
    for (int j = 0; j <= dist_x; ++j){

      // Get points for the current graticule cell, in this order:
      // bottom-left, bottom-right, top-right, top-left
      XYPoint v0;
      v0.x = cur_x;
      v0.y = cur_y;
      XYPoint v1;
      v1.x = std::min(double(lx), v0.x + 1.0);
      v1.y = v0.y;
      XYPoint v2;
      v2.x = std::min(double(lx), v0.x + 1.0);
      v2.y = std::min(double(ly), v0.y + 1.0);
      XYPoint v3;
      v3.x = v0.x;
      v3.y = std::min(double(ly), v0.y + 1.0);

      // XYPoint v0;
      // v0.x = double(int_x) + 0.5;
      // v0.y = double(int_y) + 0.5;
      // XYPoint v1;
      // v1.x = v0.x + 1.0;
      // v1.y = v0.y;
      // XYPoint v2;
      // v2.x = v0.x + 1.0;
      // v2.y = v0.y + 1.0;
      // XYPoint v3;
      // v3.x = v0.x;
      // v3.y = v0.y + 1.0;

      std::vector<XYPoint> graticule_intersections;

      // Bottom intersection
      graticule_intersections.push_back(calc_intersection(a, b, v0, v1));

      // Left intersection
      graticule_intersections.push_back(calc_intersection(a, b, v0, v3));

      // Right intersection
      graticule_intersections.push_back(calc_intersection(a, b, v1, v2));

      // Top intersection
      graticule_intersections.push_back(calc_intersection(a, b, v3, v2));

      // Diagonal intersections
      graticule_intersections.push_back(calc_intersection(a, b, v0, v2));
      graticule_intersections.push_back(calc_intersection(a, b, v3, v1));

      // Add only those intersections that are between `a` and `b`. Usually,
      // it is enough to check that the x-coordinate of the intersection is
      // between a.x and b.x. However, in some edge cases, it is possible that
      // the x-coordinate is between a.x and b.x, but the y coordinate
      // is not between a.y and b.y (e.g. if the line from a to b is
      // vertical).
      for (XYPoint inter : graticule_intersections) {
        if (((a.x <= inter.x && inter.x <= b.x) ||
             (b.x <= inter.x && inter.x <= a.x)) &&
            ((a.y <= inter.y && inter.y <= b.y) ||
             (b.y <= inter.y && inter.y <= a.y))) {
          temp_intersections.push_back(inter);
        }
      }

      cur_x += (cur_x == 0.0) ? 0.5 : 1.0;
    }
    cur_x = start_v.x;
    cur_y += (cur_y == 0.0) ? 0.5 : 1.0;
  }

  // Sort intersections
  std::sort(temp_intersections.begin(), temp_intersections.end());

  // Reverse if needed
  if ((pt1[0] > pt2[0]) || ((pt1[0] == pt2[0]) && (pt1[1] > pt2[1]))) {
    std::reverse(temp_intersections.begin(), temp_intersections.end());
  }

  // Eliminate duplicates
  std::vector<Point> intersections;
  intersections.push_back(Point(temp_intersections[0].x,
                                temp_intersections[0].y));
  for (unsigned int i = 1; i < temp_intersections.size(); ++i) {
    if (!xy_points_almost_equal(temp_intersections[i - 1],
                               temp_intersections[i])) {
      intersections.push_back(Point(temp_intersections[i].x,
                                    temp_intersections[i].y));
    }
  }
  return intersections;
}
