#include "../inset_state.h"
#include "../constants.h"
#include <CGAL/intersections.h>

// This function takes two lines as input:
// - line `a`, defined by points a1 and a2.
// - line `b`, defined by points b1 and b2.
// The function returns the intersection between them. If the two lines are
// parallel or are the same, the function returns the point (-1, -1), which
// is always outside of any graticule grid cell.
Point calc_intersection(Segment s1, Segment s2) {

  // TODO: I am unsure why this is needed. Duplicate points would be removed
  // either way.
  // Check if any segment is undefined (i.e., defined by identical points)
  // if (a1 == a2 || b1 == b2) {
  //   std::cerr << "ERROR: End points of line segment are identical"
  //             << std::endl;
  //   _Exit(EXIT_FAILURE);
  // }

  // Return segment intersection, if it exists
  const auto result = CGAL::intersection(s1, s2);
  if (result) {
    if (const Point* p = boost::get<Point>(&*result)) {
      return (*p);
    }
  }

  // Value outside grid cell indicating that no intersection was found
  return Point(-1, -1);
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
std::vector<Point> densification_points(Point pt1, Point pt2)
{
  // Vector for storing intersections before removing duplicates
  std::vector<Point> temp_intersections;

  // Assign `a` as the leftmost point of p1 and pt2. If both points have the
  // same x-coordinate, assign `a` as the lower point. The other point is
  // assigned as b. The segments (a, b) and (b, a) describe the same segment.
  // However, if we flip the order of a and b, the resulting intersections are
  // not necessarily the same because of floating point errors.
  Point a = pt1;
  Point b = pt2;
  if ((pt1[0] > pt2[0]) || ((pt1[0] == pt2[0]) && (pt1[1] > pt2[1]))) {
    std::swap(a, b);
    a = pt2;
    b = pt1;
  }
  temp_intersections.push_back(a);
  temp_intersections.push_back(b);
  Segment s1(a, b);

  // Get bottom-left point of graticule cell containing `a`
  const Point av0(floor(a.x() + 0.5) - 0.5, floor(a.y() + 0.5) - 0.5);

  // Get bottom-left point of graticule cell containing `b`
  const Point bv0(floor(b.x() + 0.5) - 0.5, floor(b.y() + 0.5) - 0.5);

  // Get bottom-left (start_v) and top-right (end_v) graticule cells of the
  // graticule cell rectangle (the smallest rectangular section of the
  // graticule grid cell containing both points)
  double y1 = av0.y();
  double y2 = bv0.y();
  if (a.y() > b.y()) {
    std::swap(y1, y2);
  }
  Point start_v(av0.x(), y1);
  Point end_v(bv0.x(), y2);

  // Loop through each row, from bottom to top
  for (unsigned int int_y = static_cast<unsigned int>(start_v.y());
       int_y <= static_cast<unsigned int>(end_v.y());
       ++int_y) {

    // Loop through each column, from left to right
    for (unsigned int int_x = static_cast<unsigned int>(start_v.x());
         int_x <= static_cast<unsigned int>(end_v.x());
         ++int_x) {

      // Get points for the current graticule cell, in this order:
      // bottom-left, bottom-right, top-right, top-left
      const Point v0(int_x + 0.5, int_y + 0.5);
      const Point v1(v0.x() + 1.0, v0.y());
      const Point v2(v0.x() + 1.0, v0.y() + 1.0);
      const Point v3(v0.x(), v0.y() + 1.0);
      std::vector<Point> graticule_intersections;

      // Bottom intersection
      graticule_intersections.push_back(calc_intersection(s1, Segment(v0, v1)));

      // Left intersection
      graticule_intersections.push_back(calc_intersection(s1, Segment(v0, v3)));

      // Right intersection
      graticule_intersections.push_back(calc_intersection(s1, Segment(v1, v2)));

      // Top intersection
      graticule_intersections.push_back(calc_intersection(s1, Segment(v3, v2)));

      // Diagonal intersections
      graticule_intersections.push_back(calc_intersection(s1, Segment(v0, v2)));
      graticule_intersections.push_back(calc_intersection(s1, Segment(v3, v1)));

      // Add only those intersections that are between `a` and `b`. Usually,
      // it is enough to check that the x-coordinate of the intersection is
      // between a.x and b.x. However, in some edge cases, it is possible that
      // the x-coordinate is between a.x and b.x, but the y coordinate
      // is not between a.y and b.y (e.g. if the line from a to b is
      // vertical).
      for (Point inter : graticule_intersections) {
        if (((a.x() <= inter.x() && inter.x() <= b.x()) ||
             (b.x() <= inter.x() && inter.x() <= a.x())) &&
            ((a.y() <= inter.y() && inter.y() <= b.y()) ||
             (b.y() <= inter.y() && inter.y() <= a.y()))) {
          temp_intersections.push_back(inter);
        }
      }
    }
  }

  // Sort intersections
  std::sort(temp_intersections.begin(), temp_intersections.end());

  // Reverse if needed
  if ((pt1[0] > pt2[0]) || ((pt1[0] == pt2[0]) && (pt1[1] > pt2[1]))) {
    std::reverse(temp_intersections.begin(), temp_intersections.end());
  }

  // Eliminate duplicates
  std::vector<Point> intersections;
  intersections.push_back(temp_intersections[0]);
  for (unsigned int i = 1; i < temp_intersections.size(); ++i) {
    if (temp_intersections[i - 1] != temp_intersections[i]) {
      intersections.push_back(temp_intersections[i]);
    }
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
      const Polygon outer = pwh.outer_boundary();
      Polygon outer_dens;

      // Iterate over each point in the outer boundary of the polygon
      for (size_t i = 0; i < outer.size(); ++i) {

        // The segment defined by points `a` and `b` is to be densified.
        // `b` should be the point immediately after `a`, unless `a` is the
        // final point of the boundary, in which case `b` should be the first
        // point.
        const Point a = outer[i];
        const Point b = (i == outer.size() - 1) ? outer[0] : outer[i + 1];

        // Densify the segment
        const std::vector<Point> outer_pts_dens = densification_points(a, b);

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
        for (size_t j = 0; j < (*h).size(); ++j) {

          // `c` and `d` are determined in the same way as `a` and `b` above
          const Point c = (*h)[j];
          const Point d = (j == (*h).size() - 1) ? (*h)[0] : (*h)[j + 1];
          const std::vector<Point> hole_pts_dens = densification_points(c, d);
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
}
