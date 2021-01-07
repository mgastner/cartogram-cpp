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

// Returns intersection of 2 lines, line a1 b1 & line a2 b2
Point calc_intersection(Point a1, Point b1, Point a2, Point b2) {

  // From https://www.geeksforgeeks.org/program-for-point-of-intersection-of-two-lines/
  // Line AB represented as a1x + b1y = c1
  double a1 = b1[1] - a1[1];
  double b1 = a1[0] - b1[0];
  double c1 = a1 * (a1[0]) + b1 * (a1[1]);

  // Line CD represented as a2x + b2y = c2
  double a2 = b2[1] - a2[1];
  double b2 = a2[0] - b2[0];
  double c2 = a2 * (a2[0]) + b2 * (a2[1]);

  double determinant = a1*b2 - a2*b1;

  if (determinant == 0)
  {
    // The lines are parallel.
    Point temp(0, 0);
    return temp;
  }
  else
  {
    double x = (b2 * c1 - b1 * c2) / determinant;
    double y = (a1 * c2 - a2 * c1) / determinant;
    Point temp(x, y);
    return temp;
  }
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

    // X of "a" < x of "b"
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

    // y of "a" < y of "b"
    for (double i = half_ceil(a[1]); i < b[1]; ++i) {
      if (i != a[1]) {
        // Get x coordinate, y coordinate = i
        double x = calc_x_intersection(a, b, i);
        Point temp(x, i);

        // Ensuring no corner points are pushed back
        if (x != half_floor(x) || (x == i && a[0] == b[0])) {
          intersections.push_back(temp);
        }
      }
    }
  } else if (a[1] > b[1]) {

    // y of "a" > y of "b"
    for (double i = half_floor(a[1]); i > b[1]; i -= 0.5) {
      if (i != a[1]) {

        // Get x coordinate, y coordinate = i
        double x = calc_x_intersection(a, b, i);
        Point temp(x, i);
        // Ensuring no corner points are pushed back
        if (x != half_floor(x) || (x == i && a[0] == b[0])) {
          intersections.push_back(temp);
        }
      }
    }
  }

  // Sorting intersections
  std::sort(intersections.begin(), intersections.end());

  // Inserting intersection with diagonals
  if (a[0] < b[0]) {
    for (unsigned int i = 0; i < intersections.size() - 1; ++i) {
      Point bottom_right(half_floor(intersections[i][0]) + 1,
                         half_floor(intersections[i][1])); // bottom right corner
      Point top_left(half_ceil(intersections[i][0]) - 1,
                     half_ceil(intersections[i][1]));; // top left corner


    }
  } else if (a[0] > b[0]) {

  } else {

  }

  // // Removing duplicates
  // intersections.erase(unique(intersections.begin(), intersections.end()),
  // intersections.end());
  if (a[0] > b[0] || (a[1] > b[1] && a[1] == b[1])) {

    // Sort in descending order of x coordinates
    std::reverse(intersections.begin(), intersections.end());
  }
  return intersections;
}
