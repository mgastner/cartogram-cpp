#include "../inset_state.h"
#include "../constants.h"
#include "../round_point.h"

// This function takes two lines as input:
// - line `a`, defined by points a1 and a2.
// - line `b`, defined by points b1 and b2.
// The function returns the intersection between them. If the two lines are
// parallel or are the same, the function returns the point (-1, -1), which
// is always outside of any graticule grid cell.
XYPoint calc_intersection(const XYPoint a1,
                          const XYPoint a2,
                          const XYPoint b1,
                          const XYPoint b2)
{
  // Check whether any segment is undefined (i.e., defined by identical
  // points)
  if (a1 == a2 || b1 == b2) {
    std::cerr << "ERROR: End points of line segment are identical"
              << std::endl;
    _Exit(EXIT_FAILURE);
  }

  // Get line equations
  const double a = (a1.y - a2.y) / (a1.x - a2.x);
  const double a_intercept = a1.y - (a1.x * a);
  const double b = (b1.y - b2.y) / (b1.x - b2.x);
  const double b_intercept = b1.y - (b1.x * b);
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

void add_intersection(std::set<XYPoint, decltype(xy_point_lesser)*>
                        *intersections,
                      const XYPoint a,
                      const XYPoint b,
                      const XYPoint c,
                      const XYPoint d,
                      const unsigned int lx,
                      const unsigned int ly)
{
  XYPoint inter = calc_intersection(a, b, c, d);
  if (((a.x <= inter.x && inter.x <= b.x) ||
        (b.x <= inter.x && inter.x <= a.x)) &&
      ((a.y <= inter.y && inter.y <= b.y) ||
        (b.y <= inter.y && inter.y <= a.y)) &&
      ((inter.y >= 0.5 && inter.y <= (ly - 0.5)) ||
        (inter.x >= 0.5 && inter.x <= (lx - 0.5)))) {
    (*intersections).insert(rounded_XYpoint(inter, lx, ly));
  }
  return;
}

void add_edge_intersection(std::set<XYPoint, decltype(xy_point_lesser)*>
                             *intersections,
                           const XYPoint a,
                           const XYPoint b,
                           const XYPoint c,
                           const XYPoint d,
                           const unsigned int lx,
                           const unsigned int ly)
{
  XYPoint inter = calc_intersection(a, b, c, d);
  if (((a.x <= inter.x && inter.x <= b.x) ||
        (b.x <= inter.x && inter.x <= a.x)) &&
      ((a.y <= inter.y && inter.y <= b.y) ||
        (b.y <= inter.y && inter.y <= a.y)) &&
      ((inter.y < 0.5 || inter.y > (ly - 0.5)) ||
        (inter.x < 0.5 || inter.x > (lx - 0.5)))) {
    (*intersections).insert(rounded_XYpoint(inter, lx, ly));
  }
  return;
}

void add_diagonals(double slope,
                   const XYPoint a,
                   const XYPoint b,
                   bool edge,
                   const unsigned int lx,
                   const unsigned int ly,
                   double start_offset,
                   std::set<XYPoint, decltype(xy_point_lesser)*>
                      *intersections)
{
  double intercept_start =
    floor(std::min(a.y - slope * a.x, b.y - slope * b.x)) + start_offset;
  double intercept_end = std::max(a.y - slope * a.x, b.y - slope * b.x);
  for (double i = intercept_start; i <= intercept_end; i += std::min(abs(slope), 1.0)) {
    XYPoint c; c.x = 0.0; c.y = i;
    XYPoint d; d.x = 1.0; d.y = slope + i;
    if (edge) {
      add_edge_intersection(intersections, a, b, c, d, lx, ly);
    } else {
      add_intersection(intersections, a, b, c, d, lx, ly);
    }
  }
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
std::vector<Point> densification_points(const Point pt1,
                                        const Point pt2,
                                        const unsigned int lx,
                                        const unsigned int ly)
{
  // If the input points are identical, return them without calculating
  // intersections
  if ((pt1.x() == pt2.x()) && (pt1.y() == pt2.y())) {
    std::vector<Point> points;
    points.push_back(pt1);
    points.push_back(pt2);
    return points;
  }

  // Ordered set for storing intersections before removing duplicates
  std::set<XYPoint, decltype(xy_point_lesser)*>
    temp_intersections(xy_point_lesser);

  // Store the leftmost point of p1 and pt2 as `a`. If both points have the
  // same x-coordinate, then store the lower point as `a`. The other point is
  // stored as `b`. The segments (a, b) and (b, a) describe the same segment.
  // However, if we flip the order of a and b, the resulting intersections are
  // not necessarily the same because of floating point errors.
  XYPoint a;
  XYPoint b;
  if ((pt1.x() > pt2.x()) || ((pt1.x() == pt2.x()) && (pt1.y() > pt2.y()))) {
    a.x = pt2.x();
    a.y = pt2.y();
    b.x = pt1.x();
    b.y = pt1.y();
  } else{
    a.x = pt1.x();
    a.y = pt1.y();
    b.x = pt2.x();
    b.y = pt2.y();
  }
  temp_intersections.insert(a);
  temp_intersections.insert(b);

  // Get vertical intersections
  double x_start = floor(a.x + 0.5) + 0.5;
  double x_end = b.x;
  for (double i = x_start; i <= x_end; i += (i == 0.0) ? 0.5 : 1.0) {
    XYPoint c; c.x = std::min(i, static_cast<double>(lx)); c.y = 0.0;
    XYPoint d; d.x = std::min(i, static_cast<double>(lx)); c.y = 1.0;
    add_intersection(&temp_intersections, a, b, c, d, lx, ly);
  }

  // Get horizontal intersections
  double y_start = floor(std::min(a.y, b.y) + 0.5) + 0.5;
  double y_end = std::max(a.y, b.y);
  for (double i = y_start; i <= y_end; i += (i == 0.0) ? 0.5 : 1.0) {
    XYPoint c; c.x = 0.0; c.y = std::min(i, static_cast<double>(ly));
    XYPoint d; d.x = 1.0; d.y = std::min(i, static_cast<double>(ly));
    add_intersection(&temp_intersections, a, b, c, d, lx, ly);
  }

  // Get bottom-left to top-right diagonal intersections
  add_diagonals(1.0, a, b, false, lx, ly, 1.0, &temp_intersections);

  // Get top-left to bottom-right diagonal intersections
  add_diagonals(-1.0, a, b, false, lx, ly, 1.0, &temp_intersections);

  // Add edge diagonals when at least one point is on the edge of the grid.
  if (a.x < 0.5 || b.x < 0.5 || a.x > (lx - 0.5) || b.x > (lx - 0.5)){

    // Bottom-left to top-right edge diagonals
    add_diagonals(2.0, a, b, true, lx, ly, 0.5, &temp_intersections);

    // Top-left to bottom-right edge diagonals
    add_diagonals(-2.0, a, b, true, lx, ly, 0.5, &temp_intersections);
  }
  if (a.y < 0.5 || b.y < 0.5 || a.y > (ly - 0.5) || b.y > (ly - 0.5)){

    // Bottom-left to top-right edge diagonals
    add_diagonals(0.5, a, b, true, lx, ly, 0.25, &temp_intersections);

    // Top-left to botom-right edge diagonals
    add_diagonals(-0.5, a, b, true, lx, ly, 0.25, &temp_intersections);
  }

  // double intercept_start = floor(std::min(a.y - a.x, b.y - b.x)) + 1.0;
  // double intercept_end = std::max(a.y - a.x, b.y - b.x);
  // for (double i = intercept_start; i <= intercept_end; ++i) {
  //   XYPoint c; c.x = 0.0; c.y = i;
  //   XYPoint d; d.x = 1.0; d.y = 1.0 + i;
  //   add_intersection(&temp_intersections, a, b, c, d, lx, ly);
  // }

  // // Get top-left to bottom-right diagonal intersections
  // intercept_start = floor(std::min(a.y + a.x, b.y + b.x)) + 1.0;
  // intercept_end = std::max(a.y + a.x, b.y + b.x);
  // for (double i = intercept_start; i <= intercept_end; ++i) {
  //   XYPoint c; c.x = 0.0; c.y = i;
  //   XYPoint d; d.x = 1.0; d.y = -1.0 + i;
  //   add_intersection(&temp_intersections, a, b, c, d, lx, ly);
  // }

  // double intercept_start;
  // double intercept_end;
  
  // // Add edge diagonals when at least one point is on the edge of the grid.
  // if (a.x < 0.5 || b.x < 0.5 || a.x > (lx - 0.5) || b.x > (lx - 0.5)){

  //   // Bottom-left to top-right edge diagonals
  //   intercept_start = floor(std::min(a.y - 2 * a.x, b.y - 2 * b.x) + 0.5) + 0.5;
  //   intercept_end = std::max(a.y - 2 * a.x, b.y - 2 * b.x);
  //   for (double i = intercept_start; i <= intercept_end; ++i) {
  //     XYPoint c; c.x = 0.0; c.y = i;
  //     XYPoint d; d.x = 1.0; d.y = 2.0 + i;
  //     add_edge_intersection(&temp_intersections, a, b, c, d, lx, ly);
  //   }

  //   // Top-left to bottom-right edge diagonals
  //   intercept_start = floor(std::min(a.y + 2 * a.x, b.y + 2 * b.x) + 0.5) + 0.5;
  //   intercept_end = std::max(a.y + 2 * a.x, b.y + 2 * b.x);
  //   for (double i = intercept_start; i <= intercept_end; ++i) {
  //     XYPoint c; c.x = 0.0; c.y = i;
  //     XYPoint d; d.x = 1.0; d.y = -2.0 + i;
  //     add_edge_intersection(&temp_intersections, a, b, c, d, lx, ly);
  //   }
  // }
  // if (a.y < 0.5 || b.y < 0.5 || a.y > (ly - 0.5) || b.y > (ly - 0.5)){

  //   // Bottom-left to top-right edge diagonals
  //   intercept_start = floor(std::min(a.y - 0.5 * a.x, b.y - 0.5 * b.x) + 0.5) - 0.25;
  //   intercept_end = std::max(a.y - 0.5 * a.x, b.y - 0.5 * b.x);
  //   for (double i = intercept_start; i <= intercept_end; i += 0.5) {
  //     XYPoint c; c.x = 0.0; c.y = i;
  //     XYPoint d; d.x = 1.0; d.y = 0.5 + i;
  //     add_edge_intersection(&temp_intersections, a, b, c, d, lx, ly);
  //   }

  //   // Top-left to botom-right edge diagonals
  //   intercept_start = floor(std::min(a.y + 0.5 * a.x, b.y + 0.5 * b.x) + 0.5) - 0.25;
  //   intercept_end = std::max(a.y + 0.5 * a.x, b.y + 0.5 * b.x);
  //   for (double i = intercept_start; i <= intercept_end; i += 0.5) {
  //     XYPoint c; c.x = 0.0; c.y = i;
  //     XYPoint d; d.x = 1.0; d.y = -0.5 + i;
  //     add_edge_intersection(&temp_intersections, a, b, c, d, lx, ly);
  //   }
  // }

  // // DEBUGGING: Check if there are any two almost equal points in the set
  // std::vector<XYPoint> inter_test(temp_intersections.begin(),
  //                                 temp_intersections.end());
  // for (unsigned int i = 1; i < inter_test.size(); ++i){
  //   if (xy_points_almost_equal(inter_test[i - 1], inter_test[i]))
  //     std::cout << "Almost equal points in set!\n";
  // }

  // Convert the set of XYPoints to a vector of CGAL points
  // TO-DO: Phase out XYPoints entirely
  std::vector<Point> intersections;
  for (auto xypt : temp_intersections)
    intersections.push_back(Point(xypt.x, xypt.y));

  // Reverse if needed
  if ((pt1.x() > pt2.x()) || ((pt1.x() == pt2.x()) && (pt1.y() > pt2.y()))) {
    std::reverse(intersections.begin(), intersections.end());
  }
  return intersections;
}

// TODO: This function may be more meaningfully included in
// check_topology.cpp.
bool duplicates(std::vector<Point> v) {
  CGAL::set_pretty_mode(std::cerr);
  for (size_t i = 0; i < v.size() - 1; ++i) {
    if (points_almost_equal(v[i], v[i + 1])) {
      std::cerr << "i = " << i << std::endl;
      std::cerr << "Point: " << i << ", v[i]: " << v[i] << std::endl;
      std::cerr << "Point: "
                << i + 1
                << ", v[i + 1]: "
                << v[i + 1]
                << std::endl;
      return true;
    }
  }
  return false;
}

void InsetState::densify_geo_divs()
{
  std::cerr << "Densifying" << std::endl;
  std::vector<GeoDiv> geodivs_dens;
  for (const auto &gd : geo_divs_) {
    GeoDiv gd_dens(gd.id());
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto outer = pwh.outer_boundary();
      Polygon outer_dens;

      // Iterate over each point in the outer boundary of the polygon
      for (size_t i = 0; i < outer.size(); ++i) {

        // The segment defined by points `a` and `b` is to be densified.
        // `b` should be the point immediately after `a`, unless `a` is the
        // final point of the boundary, in which case `b` should be the first
        // point.
        const auto a = outer[i];
        const auto b = (i == outer.size() - 1) ? outer[0] : outer[i + 1];

        // Densify the segment
        const std::vector<Point> outer_pts_dens =
          densification_points(a, b, lx_, ly_);

        // Push all points. Omit the last point because it will be included
        // in the next iteration. Otherwise, we would have duplicated points
        // in the polygon.
        for (size_t i = 0; i < (outer_pts_dens.size() - 1); ++i) {
          outer_dens.push_back(outer_pts_dens[i]);
        }
      }

      // Check for duplicate points in the densified outer boundary
      // std::vector<Point> temp_out;
      // for (size_t i = 0; i < outer_dens.size(); ++i) {
      //   temp_out.push_back(outer_dens[i]);
      // }
      // if (duplicates(temp_out)) {
      //   std::cerr << "Duplicates found in outer boundary!" << std::endl;
      // }

      std::vector<Polygon> holes_v_dens;

      // Iterate over each hole
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        Polygon hole_dens;
        for (size_t j = 0; j < h->size(); ++j) {

          // `c` and `d` are determined in the same way as `a` and `b` above
          const Point c = (*h)[j];
          const Point d = (j == h->size() - 1) ? (*h)[0] : (*h)[j + 1];
          const std::vector<Point> hole_pts_dens =
            densification_points(c, d, lx_, ly_);
          for (size_t i = 0; i < (hole_pts_dens.size() - 1); ++i) {
            hole_dens.push_back(hole_pts_dens[i]);
          }
        }

        // Check for duplicate points in the densified hole boundary
        // std::vector<Point> temp_holes;
        // for (size_t i = 0; i < hole_dens.size(); ++i) {
        //   temp_holes.push_back(hole_dens[i]);
        // }
        // if (duplicates(temp_holes)) {
        //   std::cerr << "Duplicates found in hole!" << std::endl;
        // }

        holes_v_dens.push_back(hole_dens);
      }
      const Polygon_with_holes pwh_dens(outer_dens,
                                        holes_v_dens.begin(),
                                        holes_v_dens.end());
      gd_dens.push_back(pwh_dens);
    }
    geodivs_dens.push_back(gd_dens);
  }
  set_geo_divs(geodivs_dens);
  return;
}
