#include "map_state.h"

// Returns ceiling up to nearest 0.5 value, e.g. 2.34 returns 2.5
double half_ceil(double num) {
  double decimal = num - floor(num);
  if (decimal > 0.5) {
    return ceil(num) + 0.5;
  }
  return floor(num) + 0.5;
}

// Returns floor up to nearest 0.5 value, e.g. 2.84 returns 2.5
double half_floor(double num) {
  double decimal = num - floor(num);
  if (decimal >= 0.5) {
    return floor(num) + 0.5;
  }
  return floor(num) - 0.5;
}

// Returns intersection of 2 lines, line p1p2 & line p3p4
Point calc_intersection(Point p1, Point p2, Point p3, Point p4) {

  // From https://www.geeksforgeeks.org/program-for-point-of-intersection-of-two-lines/
  // Line p1p2 represented as a1x + b1y = c1
  double a1 = p2[1] - p1[1];
  double b1 = p1[0] - p2[0];
  double c1 = a1 * (p1[0]) + b1 * (p1[1]);

  // Line p3p4 represented as a2x + b2y = c2
  double a2 = p4[1] - p3[1];
  double b2 = p3[0] - p4[0];
  double c2 = a2 * (p3[0]) + b2 * (p3[1]);

  double determinant = a1 * b2 - a2 * b1;
  double x = (b2 * c1 - b1 * c2) / determinant;
  double y = (a1 * c2 - a2 * c1) / determinant;
  Point temp(x, y);
  return temp;
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
// 1 unit apart.
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

  // Inserting intersection with diagonals
  // Checking line not through graticule line
  if (!(a[0] == b[0] && a[0] == half_floor(a[0]))) {
    std::cout << "BR-TL Run" << '\n';
    for (unsigned int i = 0; i < size - 1; ++i) {

      // Bottom-right corner
      Point bottom_right(half_floor(intersections[i][0]) + 1,
                         half_floor(intersections[i][1]));

      // Top-left corner
      Point top_left(half_ceil(intersections[i + 1][0]) - 1,
                     half_ceil(intersections[i + 1][1]));

      // Finding intersection
      Point temp = calc_intersection(intersections[i], intersections[i + 1],
                                     top_left, bottom_right);
      intersections.push_back(temp);
    }
  }

  // Checking line not through middle of graticule line
  if (a[0] == b[0] && a[0] != floor(a[0])) {
    std::cout << "BL-TR Run" << '\n';
    for (unsigned int i = 0; i < size - 1; ++i) {

      // Bottom-left corner
      Point bottom_left(half_floor(intersections[i][0]),
                         half_floor(intersections[i][1]));

      // Top-right corner
      Point top_right(half_ceil(intersections[i + 1][0]),
                     half_ceil(intersections[i + 1][1]));

      // Finding intersection
      Point temp = calc_intersection(intersections[i], intersections[i + 1],
                                     top_right, bottom_left);
      intersections.push_back(temp);
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
