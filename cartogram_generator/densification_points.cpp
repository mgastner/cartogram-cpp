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

  // Vertical intersections (intersections with x coordinate lines)
  if (a[0] < b[0]) {

    // X of "a" < x of "b", "x" is each vertical graticule line
    for (double x = half_ceil(a[0]); x < b[0]; ++x) {
      if (x != a[0]) {

        // Get y coordinate, x coordinate = x
        double y = calc_y_intersection(a, b, x);
        Point temp(x, y);
        intersections.push_back(temp);
      }
    }
  } else if (a[0] > b[0]) {

    // x of "a" > x of "b"
    for (double x = half_floor(a[0]); x > b[0]; --x) {
      if (x != a[0]) {

        // Get y coordinate, x coordinate = x
        double y = calc_y_intersection(a, b, x);
        Point temp(x, y);
        intersections.push_back(temp);
      }
    }
  }

  // Horizontal intersections (intersections with y coordinate lines)
  if (a[1] < b[1]) {

    // Y of "a" < y of "b", "y" is each horizontal graticule line
    for (double y = half_ceil(a[1]); y < b[1]; ++y) {
      if (y != a[1]) {
        // Get x coordinate, y coordinate = y
        double x = calc_x_intersection(a, b, y);
        Point temp(x, y);

        // Ensuring no corner points are pushed back
        // As they would've already been added when calculating
        // Vertical intersections
        if (x != half_floor(x) || a[0] == b[0]) {
          intersections.push_back(temp);
        }
      }
    }
  } else if (a[1] > b[1]) {

    // y of "a" > y of "b"
    for (double y = half_floor(a[1]); y > b[1]; --y) {
      if (y != a[1]) {

        // Get x coordinate, y coordinate = y
        double x = calc_x_intersection(a, b, y);
        Point temp(x, y);
        // Ensuring no corner points are pushed back
        if (x != half_floor(x) || a[0] == b[0]) {
          intersections.push_back(temp);
        }
      }
    }
  }

  // std::cout << "Graticule intersections calculated; Now sorting..." << '\n';

  // Sorting intersections
  std::sort(intersections.begin(), intersections.end());

  // std::cout << "Sorted!\n";;

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
    Point top_left(half_ceil(intersections[i + 1][0]) - 1,
                   half_ceil(intersections[i + 1][1]));

    // Finding intersection

    // std::cout << "Calculating segments..." << '\n';

    Segment seg_diaganol(top_left, bottom_right);
    Segment seg_intersec(intersections[i], intersections[i + 1]);

    auto result = intersection(seg_diaganol, seg_intersec);
    if (result) {
      if (const Point *s = boost::get<Point>(&*result)) {
        intersections.push_back(*s);
      }
    }
  }

  // Adding intersections with bottom_left-top_right diaganol
  for (unsigned int i = 0; i < size - 1; ++i) {

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

        // Checking intersection not through middle of graticule lines,
        // Because if middle of graticule lines, then duplicates will be added as
        // intersection already added by previous for loop.
        if (!((*s)[0] == floor((*s)[0]) || (*s)[1] == floor((*s)[1]))) {
          intersections.push_back(*s);
        }
      }
    }
  }

  if (size == 0) {

    // Adding intersection with singular grid cell diaganols

    // Bottom-right corner
    Point bottom_right(half_floor(a[0]) + 1,
                       half_floor(a[1]));

    // Top-left corner
    Point top_left(half_ceil(a[0]) - 1,
                   half_ceil(a[1]));

    // Bottom-left corner
    Point bottom_left(half_floor(a[0]),
                      half_floor(a[1]));

    // Top-right corner
    Point top_right(half_ceil(a[0]),
                    half_ceil(a[1]));

    // Finding intersections
    Segment seg_tlbr(top_left, bottom_right);
    Segment seg_trbl(top_right, bottom_left);
    Segment seg_intersec(a, b);

    auto result_1 = intersection(seg_tlbr, seg_intersec);
    if (result_1) {
      if (const Point *s = boost::get<Point>(&*result_1)) {
        intersections.push_back(*s);
      }
    }

    auto result_2 = intersection(seg_trbl, seg_intersec);
    if (result_2) {
      if (const Point *s = boost::get<Point>(&*result_2)) {
        if (!((*s)[0] == floor((*s)[0]) || (*s)[1] == floor((*s)[1]))) {
          intersections.push_back(*s);
        }
      }
    }

  } else {

    // Adding final intersection using last intersection point

    // Bottom-right corner
    Point bottom_right(half_floor(a[0]) + 1,
                       half_floor(a[1]));

    // Top-left corner
    Point top_left(half_ceil(a[0]) - 1,
                   half_ceil(a[1]));

    // Bottom-left corner
    Point bottom_left(half_floor(a[0]),
                      half_floor(a[1]));

    // Top-right corner
    Point top_right(half_ceil(a[0]),
                    half_ceil(a[1]));

    // Finding intersections
    Segment seg_tlbr(top_left, bottom_right);
    Segment seg_trbl(top_right, bottom_left);

    Point c;
    if (a[0] > b[0] || (a[1] > b[1] && a[0] == b[0])) {
      c = a;
    } else {
      c = b;
    }

    Segment seg_intersec(intersections[size - 1], c);

    auto result_1 = intersection(seg_tlbr, seg_intersec);
    if (result_1) {
      if (const Point *s = boost::get<Point>(&*result_1)) {
        intersections.push_back(*s);
      }
    }

    auto result_2 = intersection(seg_trbl, seg_intersec);
    if (result_2) {
      if (const Point *s = boost::get<Point>(&*result_2)) {
        if (!((*s)[0] == floor((*s)[0]) || (*s)[1] == floor((*s)[1]))) {
          intersections.push_back(*s);
        }
      }
    }
  }

  // Removing duplicates
  // intersections.erase(unique(intersections.begin(), intersections.end()),
  // intersections.end());

  // std::cout << "Diaganol intersections calculated! Now sorting..." << '\n';

  // Sorting intersections
  std::sort(intersections.begin(), intersections.end());
  if (a[0] > b[0] || (a[1] > b[1] && a[0] == b[0])) {

    // Sort in descending order of x coordinates
    std::reverse(intersections.begin(), intersections.end());
  }

  // std::cout << "Sorted! Now returning..." << '\n';

  return intersections;
}
