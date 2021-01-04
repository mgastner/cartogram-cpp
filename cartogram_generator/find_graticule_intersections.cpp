#include "map_state.h"

// This function takes two points, "a" and "b", and returns all horizontal and
// vertical intersections with a graticule with graticule lines placed
// 1 unit apart.
std::vector<Point> graticule_intersections(Point a, Point b)
{
  std::vector<Point> intersections;

  // Vertical intersections (x coordinate lines)
  if (a[0] < b[0]) {

    // x of "a" < x of "b"
    for (unsigned int i = ceil(a[0]); i < b[0]; ++i) {
      if (i != a[0]) {

        // get y coordinate, x coordinate = i
        double y = (a[1] * (b[0] - i) +
                    b[1] * (i - a[0])) /
                   (b[0] - a[0]);
        Point temp(i, y);
        intersections.push_back(temp);
      }
    }
  } else if (a[0] > b[0]) {

    // x of "a" > x of "b"
    for (unsigned int i = floor(a[0]); i > b[0]; --i) {
      if (i != a[0]) {

        // get y coordinate, x coordinate = i
        double y = (a[1] * (b[0] - i) +
                    b[1] * (i - a[0])) /
                   (b[0] - a[0]);
        Point temp(i, y);
        intersections.push_back(temp);
      }
    }
  }

  // Horizontal intersections (y coordinate lines)
  if (a[1] < b[1]) {

    // y of "a" < y of "b"
    for (unsigned int i = ceil(a[1]); i < b[1]; ++i) {
      if (i != a[1]) {
        // get x coordinate, y coordinate = i
        double x = (a[0] * (b[1] - i) +
                    b[0] * (i - a[1])) /
                   (b[1] - a[1]);
        Point temp(x, i);

        // Ensuring no corner points are pushed back
        if (x != i || (x == i && a[0] == b[0])) {
          intersections.push_back(temp);
        }
      }
    }
  } else if (a[1] > b[1]) {

    // y of "a" > y of "b"
    for (unsigned int i = floor(a[1]); i > b[1]; --i) {
      if (i != a[1]) {

        // get x coordinate, y coordinate = i
        double x = (a[0] * (b[1] - i) +
                    b[0] * (i - a[1])) /
                   (b[1] - a[1]);
        Point temp(x, i);
        // Ensuring no corner points are pushed back
        if (x != i || (x == i && a[0] == b[0])) {
          intersections.push_back(temp);
        }
      }
    }
  }

  // Sorting intersections
  std::sort(intersections.begin(), intersections.end());
  if (a[0] > b[0] || (a[1] > b[1] && a[1] == b[1])) {

    // Sort in descending order of x coordinates
    std::reverse(intersections.begin(), intersections.end());
  }
  return intersections;
}
