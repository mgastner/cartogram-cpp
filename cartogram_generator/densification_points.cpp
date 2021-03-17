#include "map_state.h"
#include <CGAL/intersections.h>

// Returns ceiling up to nearest 0.5 value, e.g. 2.64 returns 3.5
double half_ceil(double num) {
  return ceil(num - 0.5) + 0.5;
}

// Returns floor up to nearest 0.5 value, e.g. 2.23 returns 1.5
double half_floor(double num) {
  return floor(num + 0.5) - 0.5;
}

// Returns intersection of line with vertical grid line
double calc_y_intersection(Point a, Point b, double i) {
  return (a[1] * (b[0] - i) + b[1] * (i - a[0])) / (b[0] - a[0]);
}

// Returns intersection of line with horizontal grid line
double calc_x_intersection(Point a, Point b, double i) {
  return (a[0] * (b[1] - i) + b[0] * (i - a[1])) / (b[1] - a[1]);
}

// This function takes two points, "a" and "b", and returns all horizontal and
// vertical intersections with a graticule with graticule lines placed
// 1 unit apart. Further, it returns all points that collide with the diaganol
// of grid cells. Assumes grid cells to start at 0.5, 0.5
std::vector<Point> densification_points(Point a, Point b)
{
  std::vector<Point> intersections;

  // Vertical intersections (x coordinate lines)
  if (a[0] < b[0]) {

    // X of "a" < x of "b", "i" is each vertical graticule line
    for (double i = half_ceil(a[0]); i < b[0]; ++i) {
      if (i != a[0]) {

        // Get y coordinate, x coordinate = i
        double y = calc_y_intersection(a, b, i);
        Point temp(i, y);
        intersections.push_back(temp);
      }
    }
  } else if (a[0] > b[0]) {

    // x of "a" > x of "b"
    for (double i = half_floor(a[0]); i > b[0]; --i) {
      if (i != a[0]) {

        // Get y coordinate, x coordinate = i
        double y = calc_y_intersection(a, b, i);
        Point temp(i, y);
        intersections.push_back(temp);
      }
    }
  }

  // Horizontal intersections (y coordinate lines)
  if (a[1] < b[1]) {

    // y of "a" < y of "b", "i" is each horizontal graticule line
    for (double i = half_ceil(a[1]); i < b[1]; ++i) {
      if (i != a[1]) {
        // Get x coordinate, y coordinate = i
        double x = calc_x_intersection(a, b, i);
        Point temp(x, i);

        // Ensuring no corner points are pushed back
        if (x != half_floor(x) || a[0] == b[0]) {
          intersections.push_back(temp);
        }
      }
    }
  } else if (a[1] > b[1]) {

    // y of "a" > y of "b"
    for (double i = half_floor(a[1]); i > b[1]; --i) {
      if (i != a[1]) {

        // Get x coordinate, y coordinate = i
        double x = calc_x_intersection(a, b, i);
        Point temp(x, i);
        // Ensuring no corner points are pushed back
        if (x != half_floor(x) || a[0] == b[0]) {
          intersections.push_back(temp);
        }
      }
    }
  }

  // Sorting intersections
  std::sort(intersections.begin(), intersections.end());

  // Storing current size
  unsigned int size = intersections.size();

  // TODO: Fix bug for diaganols when 0 intersections

  // Inserting intersection with diagonals
  for (unsigned int i = 0; i < size + 1; ++i) {

    // Bottom-right corner
    Point bottom_right(half_floor(intersections[i][0]) + 1,
                       half_floor(intersections[i][1]));

    // Top-left corner
    Point top_left(half_ceil(intersections[i + 1][0]) - 1,
                   half_ceil(intersections[i + 1][1]));

    // Finding intersection

    Segment seg_diaganol(top_left, bottom_right);
    Segment seg_intersec(intersections[i], intersections[i + 1]);

    auto result = intersection(seg_diaganol, seg_intersec);
    if (result) {
      if (const Point *s = boost::get<Point>(&*result)) {
        intersections.push_back(*s);
      }
    }
  }

  // Checking line not through middle of graticule lines,
  // Because if middle of graticule lines, then duplicates will be added as
  // intersection already accounted for by previous for loop.
  if (!(a[0] == b[0] && a[0] == floor(a[0])) &&
      !(a[1] == b[1] && a[1] == floor(a[1]))) {
    for (unsigned int i = 0; i < size + 1; ++i) {

      // Bottom-left corner
      Point bottom_left(half_floor(intersections[i][0]),
                         half_floor(intersections[i][1]));

      // Top-right corner
      Point top_right(half_ceil(intersections[i + 1][0]),
                     half_ceil(intersections[i + 1][1]));

      // Finding intersection
      Segment seg_diaganol(top_right, bottom_left);
      Segment seg_intersec(intersections[i], intersections[i + 1]);

      auto result = intersection(seg_diaganol, seg_intersec);
      if (result) {
        if (const Point *s = boost::get<Point>(&*result)) {
          intersections.push_back(*s);
        }
      }
    }
  }

  // Removing duplicates
  // intersections.erase(unique(intersections.begin(), intersections.end()),
  // intersections.end());

  // Sorting intersections
  std::sort(intersections.begin(), intersections.end());
  if (a[0] > b[0] || (a[1] > b[1] && a[1] == b[1])) {

    // Sort in descending order of x coordinates
    std::reverse(intersections.begin(), intersections.end());
  }
  return intersections;
}
