#include "../inset_state.h"
#include "../constants.h"
#include <CGAL/intersections.h>
#include "../round_point.h"

// This function takes two lines as input:
// - line `a`, defined by points a1 and a2.
// - line `b`, defined by points b1 and b2.
// The function returns the intersection between them. If the two lines are
// parallel or are the same, the function returns the point (-1, -1), which
// is always outside of any graticule grid cell.
Point calc_intersection(const Point a,
                        const Point b,
                        const double coef_x,
                        const double coef_y,
                        const double intercept,
                        const unsigned int lx,
                        const unsigned int ly)
{
  const auto result = CGAL::intersection(
    Line(coef_x, coef_y, intercept),
    Segment(a, b)
  );
  if (result) {
    if (const Point* p = boost::get<Point>(&*result)) {
      return rounded_point((*p), lx, ly);
    }
  }
  return Point(-1.0, -1.0);
}

// edge: 'x', 'y', 'n'
void add_diag_intersection(std::set<Point, decltype(point_lesser)*>
                              *intersections,
                            const Point a,
                            const Point b,
                            double slope,
                            double base_intercept,
                            double step,
                            const unsigned int lx,
                            const unsigned int ly,
                            char edge)
{
  double intercept_start =
    floor(std::min(a.y() - slope * a.x(), b.y() - slope * b.x()))
    + base_intercept;
  double intercept_end = std::max(a.y() - slope * a.x(), b.y() - slope * b.x());
  for (double i = intercept_start; i <= intercept_end; i += step) {
    Point inter = calc_intersection(a, b, slope, -1.0, i, lx, ly);
    if (inter > Point(0, 0)){
      if ((edge == 'x' && (inter.x() < 0.5 || inter.x() > (lx - 0.5))) ||
          (edge == 'y' && (inter.y() < 0.5 || inter.y() > (ly - 0.5))) ||
          (edge == 'n' && inter.x() >= 0.5 && inter.x() <= (lx - 0.5) &&
            inter.y() >= 0.5 && inter.y() <= (ly - 0.5))){
        (*intersections).insert(inter);
      }
    }
  }
}

// TODO: If a or b are themselves intersection points (e.g. if pt1.x is an
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
  std::set<Point, decltype(point_lesser)*>
    temp_intersections(point_lesser);

  // Store the leftmost point of p1 and pt2 as `a`. If both points have the
  // same x-coordinate, then store the lower point as `a`. The other point is
  // stored as `b`. The segments (a, b) and (b, a) describe the same segment.
  // However, if we flip the order of a and b, the resulting intersections are
  // not necessarily the same because of floating point errors.
  Point a;
  Point b;
  if ((pt1.x() > pt2.x()) || ((pt1.x() == pt2.x()) && (pt1.y() > pt2.y()))) {
    a = pt2; b = pt1;
  } else{
    a = pt1; b = pt2;
  }
  temp_intersections.insert(a);
  temp_intersections.insert(b);

  // Get vertical intersections
  double x_start = floor(a.x() + 0.5) + 0.5;
  double x_end = b.x();
  for (double i = x_start; i <= x_end; i += (i == 0.0) ? 0.5 : 1.0) {
    Point inter = calc_intersection(a, b, 1.0, 0.0, -i, lx, ly);
    if (inter > Point(0, 0)) temp_intersections.insert(inter);
  }

  // Get horizontal intersections
  double y_start = floor(std::min(a.y(), b.y()) + 0.5) + 0.5;
  double y_end = std::max(a.y(), b.y());
  for (double i = y_start; i <= y_end; i += (i == 0.0) ? 0.5 : 1.0) {
    Point inter = calc_intersection(a, b, 0.0, 1.0, -i, lx, ly);
    if (inter > Point(0, 0)) temp_intersections.insert(inter);
  }

  // Get bottom-left to top-right diagonal intersections
  add_diag_intersection(&temp_intersections, a, b, 1.0, 0.0, 1.0, lx, ly, 'n');

  // Get top-left to bottom-right diagonal intersections
  add_diag_intersection(&temp_intersections, a, b, -1.0, 0.0, 1.0, lx, ly, 'n');

  // Add edge diagonals when at least one point is on the edge of the grid.
  if (a.x() < 0.5 || b.x() < 0.5 || a.x() > (lx - 0.5) || b.x() > (lx - 0.5)){

    // Bottom-left to top-right edge diagonals
    add_diag_intersection(&temp_intersections, a, b, 2.0, 0.5, 1.0, lx, ly, 'x');

    // Top-left to bottom-right edge diagonals
    add_diag_intersection(&temp_intersections, a, b, -2.0, 0.5, 1.0, lx, ly, 'x');
  }
  if (a.y() < 0.5 || b.y() < 0.5 || a.y() > (ly - 0.5) || b.y() > (ly - 0.5)){

    // Bottom-left to top-right edge diagonals
    add_diag_intersection(&temp_intersections, a, b, 0.5, 0.25, 0.5, lx, ly, 'y');

    // Top-left to botom-right edge diagonals
    add_diag_intersection(&temp_intersections, a, b, -0.5, 0.25, 0.5, lx, ly, 'y');
  }

  // // DEBUGGING: Check if there are any two almost equal points in the set
  // std::vector<XYPoint> inter_test(temp_intersections.begin(),
  //                                 temp_intersections.end());
  // for (unsigned int i = 1; i < inter_test.size(); ++i){
  //   if (xy_points_almost_equal(inter_test[i - 1], inter_test[i]))
  //     std::cout << "Almost equal points in set!\n";
  // }

  // Create a Point vector from the set
  std::vector<Point> intersections(temp_intersections.begin(),
                                   temp_intersections.end());
  
  // Reverse if needed
  if ((pt1.x() > pt2.x()) || ((pt1.x() == pt2.x()) && (pt1.y() > pt2.y()))) {
    std::reverse(intersections.begin(), intersections.end());
  }
  return intersections;
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
