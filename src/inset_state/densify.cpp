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

  // Vector for storing intersections before removing duplicates
  std::vector<XYPoint> temp_intersections;

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
  temp_intersections.push_back(a);
  temp_intersections.push_back(b);

  // Get bottom-left point of graticule cell containing `a`
  XYPoint av0;
  av0.x = std::max(0.0, floor(a.x + 0.5) - 0.5);
  av0.y = std::max(0.0, floor(a.y + 0.5) - 0.5);

  // Get bottom-left point of graticule cell containing `b`
  XYPoint bv0;
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

  // Distance between left-most and right-most graticule cell
  const unsigned int dist_x = std::ceil(end_v.x - start_v.x);

  // Distance between top and bottom graticule cell
  const unsigned int dist_y = std::ceil(end_v.y - start_v.y);

  // Iterator variables for tracking current graticule cell in next for-loop
  double current_graticule_x = start_v.x;
  double current_graticule_y = start_v.y;

  // TODO: IT IS INEFFICIENT TO RUN THE INTERSECTIONS OF THE LINE FROM a TO b
  // WITH THE HORIZONTAL GRID LINES AND VERTICAL GRID LINES IN EACH ITERATION
  // OF THE NESTED LOOP BELOW. IT WOULD BE BETTER TO HAVE A NON-NESTED LOOP
  // OVER EACH HORIZONTAL LINE IN THE RANGE, THEN A NON-NESTED LOOP OVER EACH
  // VERTICAL LINE. FOR THE DIAGONALS, THIS PROCEDURE IS A LITTLE BIT
  // TRICKIER; AT THE EDGES THE DIAGONALS ARE NOT STRAIGHT CONTINUATIONS OF
  // THE ADJOINING DIAGONALS. STILL, THERE IS A PROBABLY A WAY TO GET
  // INTERSECTIONS WITH THE 'MAIN' DIAGONALS ADN TREAT THE EDGE CASES
  // SEPARATELY.

  // Iterate over each row, from bottom to top
  for (unsigned int i = 0; i <= dist_y; ++i) {

    // Iterate over each column, from left to right
    for (unsigned int j = 0; j <= dist_x; ++j) {

      // Get points for the current graticule cell, in the following order:
      // bottom-left, bottom-right, top-right, top-left
      const XYPoint v0(current_graticule_x, current_graticule_y);
      const XYPoint v1(
        v0.x == 0.0 ? 0.5 : std::min(double(lx), v0.x + 1.0),
        v0.y);
      const XYPoint v2(
        v1.x,
        v0.y == 0.0 ? 0.5 : std::min(double(ly), v0.y + 1.0));
      const XYPoint v3(v0.x, v2.y);

      // Store intersections of line segment from `a` to `b` with graticule
      // lines and diagonals
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
      for (const auto &inter : graticule_intersections) {
        if (((a.x <= inter.x && inter.x <= b.x) ||
             (b.x <= inter.x && inter.x <= a.x)) &&
            ((a.y <= inter.y && inter.y <= b.y) ||
             (b.y <= inter.y && inter.y <= a.y))) {
          temp_intersections.push_back(rounded_XYpoint(inter, lx, ly));
        }
      }

      // If the current graticule cell touches the left edge, add 0.5 to
      // obtain the next graticule cell. Otherwise, add 1.0.
      current_graticule_x += (current_graticule_x == 0.0) ? 0.5 : 1.0;
    }
    current_graticule_x = start_v.x;

    // If the current row touches the bottom edge, add 0.5 to
    // obtain the next row. Otherwise, add 1.0.
    current_graticule_y += (current_graticule_y == 0.0) ? 0.5 : 1.0;
  }

  // Sort intersections
  std::sort(temp_intersections.begin(), temp_intersections.end());

  // TODO: IF temp_intersections WERE AN ORDERED SET, THERE WOULD BE NO NEED
  // TO REMOVE DUPLICATES
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
