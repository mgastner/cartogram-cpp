#include "inset_state.h"
#include "constants.h"
#include <CGAL/intersections.h>

// Returns ceiling up to nearest 0.5 value, e.g. 2.64 returns 3.5
double half_ceil(double num) {

  // Extracting decimal
  double decimal = num - floor(num);

  // Checking whether decimal is 0.5
  if (abs(decimal - 0.5) <= 1e-10) {
    decimal = 0.5;
  }

  // Checking whether to ceil
  if (decimal > 0.5) {
    return ceil(num) + 0.5;
  }

  // Else flooring
  return floor(num) + 0.5;
}

// Returns floor up to nearest 0.5 value, e.g. 2.23 returns 1.5
double half_floor(double num){
  return floor(num + 0.5) - 0.5;
}

// Very similar doubles
bool almost_equal(double a, double b) {
  return abs(a - b) <= 1e-10;
}

// Very similar points
bool point_almost_equal(Point a, Point b) {
  return (almost_equal(a[0], b[0]) && almost_equal(a[1], b[1]));
}

// Truncate a double down to round_digits digits after decimal point.
// round_digits is defined in constants.h
double round_down(double value) {
  return floor(value * round_digits) / round_digits;
}

// Truncate the coordinates of a point
Point round_point(Point p1) {
  double x1 = round_down(p1[0]);
  double y1 = round_down(p1[1]);
  Point rounded_point(x1, y1);
  return rounded_point;
}

// Returns intersection of line with vertical grid line
double calc_y_intersection(Point a, Point b, double x) {
  return (a[1] * (b[0] - x) + b[1] * (x - a[0])) / (b[0] - a[0]);
}

// Returns intersection of line with horizontal grid line
double calc_x_intersection(Point a, Point b, double y) {
  return (a[0] * (b[1] - y) + b[0] * (y - a[1])) / (b[1] - a[1]);
}

// This function takes two points, "a" and "b", and returns all horizontal and
// vertical intersections with a graticule with graticule lines placed
// 1 unit apart. Further, it returns all points that collide with the diaganol
// of grid cells. Assumes grid cells to start at 0.5, 0.5
std::vector<Point> densification_points(Point a, Point b)
{

  // Vector to store all intersections
  std::vector<Point> intersections;
  intersections.push_back(a);
  intersections.push_back(b);

  // Vertical intersections (intersections with x coordinate lines)
  if (a[0] < b[0]) {

    // X of "a" < x of "b", "x" is each vertical graticule line
    for (double x = half_ceil(a[0]); x < b[0]; ++x) {
      if (!almost_equal(x, a[0]) && !almost_equal(x, b[0])) {

        // Get y coordinate, x coordinate = x
        double y = calc_y_intersection(a, b, x);
        Point temp(x, y);
        if (!point_almost_equal(a, temp) && !point_almost_equal(b, temp)) {
          intersections.push_back(round_point(temp));
        }
      }
    }
  } else if (a[0] > b[0]) {

    // x of "a" > x of "b"
    for (double x = half_floor(a[0]); x > b[0]; --x) {
      if (!almost_equal(x, a[0]) && !almost_equal(x, b[0])) {

        // Get y coordinate, x coordinate = x
        double y = calc_y_intersection(a, b, x);
        Point temp(x, y);
        if (!point_almost_equal(a, temp) && !point_almost_equal(b, temp)) {
          intersections.push_back(round_point(temp));
        }
      }
    }
  }

  // Horizontal intersections (intersections with y coordinate lines)
  if (a[1] < b[1]) {

    // Y of "a" < y of "b", "y" is each horizontal graticule line
    for (double y = half_ceil(a[1]); y < b[1]; ++y) {
      if (!almost_equal(y, a[1]) && !almost_equal(y, b[1])) {
        // Get x coordinate, y coordinate = y
        double x = calc_x_intersection(a, b, y);
        Point temp(x, y);

        // Ensuring no corner points are pushed back as they would've already
        // been added when calculating vertical intersections
        if (x != half_floor(x) || a[0] == b[0]) {
          if (!point_almost_equal(a, temp) && !point_almost_equal(b, temp)) {
            intersections.push_back(round_point(temp));
          }
        }
      }
    }
  } else if (a[1] > b[1]) {

    // y of "a" > y of "b"
    for (double y = half_floor(a[1]); y > b[1]; --y) {
      if (!almost_equal(y, a[1]) && !almost_equal(y, b[1])) {

        // Get x coordinate, y coordinate = y
        double x = calc_x_intersection(a, b, y);
        Point temp(x, y);
        if (!point_almost_equal(a, temp) && !point_almost_equal(b, temp)) {
          intersections.push_back(round_point(temp));
        }
      }
    }
  }

  // Sorting intersections
  std::sort(intersections.begin(), intersections.end());

  // Clearing duplicates
  intersections.erase(unique(intersections.begin(), intersections.end()),
  intersections.end());

  // Storing current size
  double size = intersections.size();

  // TODO: Fix bug for diaganols when 0 intersections

  // Inserting intersection with diagonals

  // Adding intersections with bottom_right-top_left diaganol
  for (unsigned int i = 0; i < size - 1; ++i) {

    // std::cout << "Inside for loop! Calculating diaganol intersections...\n";;

    // Bottom-right corner
    Point bottom_right(half_floor(intersections[i][0]) + 1,
                       half_floor(intersections[i][1]));

    // Top-left corner
    Point top_left(bottom_right[0] - 1,
                   bottom_right[1] + 1);

    // Finding intersection
    Segment seg_diaganol(top_left, bottom_right);
    Segment seg_intersec(a, b);
    auto result = CGAL::intersection(seg_diaganol, seg_intersec);
    if (result) {
      if (const Point *s = boost::get<Point>(&*result)) {
        if (!point_almost_equal(a, (*s)) && !point_almost_equal(b, (*s))) {
          intersections.push_back(round_point((*s)));
        }
      }
    }
  }

  // Adding intersections with bottom_left-top_right diaganol
  for (unsigned int i = 0; i < size - 1; ++i) {

    // Bottom-left corner
    Point bottom_left(half_floor(intersections[i][0]),
                      half_floor(intersections[i][1]));

    // Top-right corner
    Point top_right(bottom_left[0] + 1,
                    bottom_left[1] + 1);

    // Finding intersection
    Segment seg_diaganol(top_right, bottom_left);
    Segment seg_intersec(a, b);

    auto result = CGAL::intersection(seg_diaganol, seg_intersec);
    if (result) {
      if (const Point *s = boost::get<Point>(&*result)) {

        // Checking intersection not through middle of graticule lines,
        // Because if middle of graticule lines, then duplicates will be added as
        // intersection already added by previous for loop.
        if (!almost_equal((*s)[0], floor((*s)[0]))) {
          if (!point_almost_equal(a, (*s)) && !point_almost_equal(b, (*s))) {
            intersections.push_back(round_point((*s)));
          }
        }
      }
    }
  }

  // TODO: REMOVE DUPLICATES UNTIL NO MORE

  // Sorting intersections
  std::sort(intersections.begin(), intersections.end());

  // Removing duplicates
  intersections.erase(unique(intersections.begin(), intersections.end()),
  intersections.end());

  if (a[0] > b[0] || (a[1] > b[1] && almost_equal(a[0], b[0]))) {

    // Sort in descending order of x coordinates
    std::reverse(intersections.begin(), intersections.end());
  }
  
  return intersections;
}
