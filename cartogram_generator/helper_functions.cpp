#include "map_state.h"

// Struct to store based on x value
struct x_sort {

  double x;
  double y;

  // Overload "<" operator for this data type. Idea from
  // https://stackoverflow.com/questions/4892680/sorting-a-vector-of-structs
  const bool operator < (const x_sort &rhs) const
  {
    return (x < rhs.x);
  }
  const bool operator == (const x_sort &rhs) const
  {
    return (x == rhs.x && y == rhs.y);
  }

};

// Struct to sort based on y value
struct y_sort {

  double x;
  double y;

  // Overload "<" operator for this data type. Idea from
  // https://stackoverflow.com/questions/4892680/sorting-a-vector-of-structs
  const bool operator < (const y_sort &rhs) const
  {
    return (y < rhs.y);
  }
  const bool operator == (const y_sort &rhs) const
  {
    return (x == rhs.x && y == rhs.y);
  }

};

// This function takes two points, "a" and "b", and returns all horizontal and
// vertical intersections with a graticule with graticule lines placed
// 1 unit apart.
std::vector<Point> graticule_intersections(Point a, Point b)
{

  std::vector<Point> intersections;
  std::vector<Point> intersections_sorted;

  // Vertical intersections (x coordinate lines)
  if (a[0] < b[0]) {

    // x of "a" < x of "b"
    for (unsigned int i = ceil(a[0]); i < b[0]; ++i) {
      if (i == a[0]) {
        continue;
      }

      // get y coordinate, x coordinate = i
      double y = (a[1] * (b[0] - i) +
                b[1] * (i - a[0])) /
               (b[0] - a[0]);
      Point temp(i, y);
      intersections.push_back(temp);
    }
  } else if (a[0] > b[0]) {

    // x of "a" > x of "b"
    for (unsigned int i = floor(a[0]); i > b[0]; --i) {
      if (i == a[0]) {
        continue;
      }

      // get y coordinate, x coordinate = i
      double y = (a[1] * (b[0] - i) +
                b[1] * (i - a[0])) /
               (b[0] - a[0]);
      Point temp(i, y);
      intersections.push_back(temp);
    }
  }

  // Horizontal intersections (y coordinate lines)
  if (a[1] < b[1]) {

    // y of "a" < y of "b"
    for (unsigned int i = ceil(a[1]); i < b[1]; ++i) {
      if (i == a[1]) {
        continue;
      }

      // get x coordinate, y coordinate = i
      double x = (a[0] * (b[1] - i) +
                b[0] * (i - a[1])) /
               (b[1] - a[1]);
      Point temp(x, i);
      intersections.push_back(temp);
    }
  } else if (a[1] > b[1]) {

    // y of "a" > y of "b"
    for (unsigned int i = floor(a[1]); i > b[1]; --i) {
      if (i == a[1]) {
        continue;
      }

      // get x coordinate, y coordinate = i
      double x = (a[0] * (b[1] - i) +
                b[0] * (i - a[1])) /
               (b[1] - a[1]);
      Point temp(x, i);
      intersections.push_back(temp);
    }
  }
  else {

    // Error: same point repeated
  }

  if (a[0] != b[0]) {

    // sorting by x coordinate
    std::vector<x_sort> intersections_temp;
    for (Point c : intersections) {
      x_sort temp;
      temp.x = c[0];
      temp.y = c[1];
      intersections_temp.push_back(temp);
    }

    // Sort in ascending order of x coordinates
    sort(intersections_temp.begin(), intersections_temp.end());
    intersections_temp.erase(unique(intersections_temp.begin(),
    intersections_temp.end()), intersections_temp.end());
    if (a[0] > b[0]) {

      // Sort in descending order of x coordinates
      reverse(intersections_temp.begin(), intersections_temp.end());
    }

    // Converting to sorted vector of points
    for (x_sort point : intersections_temp) {
      intersections_sorted.push_back(Point(point.x, point.y));
    }
  }
  else
  {
    // sorting by y coordinate
    std::vector<y_sort> intersections_temp;
    for (Point c : intersections) {
      y_sort temp;
      temp.x = c[0];
      temp.y = c[1];
      intersections_temp.push_back(temp);
    }

    // Sort in ascending order of y coordinates
    sort(intersections_temp.begin(), intersections_temp.end());
    intersections_temp.erase(unique(intersections_temp.begin(),
    intersections_temp.end()), intersections_temp.end());

    if (a[1] > b[1]) {

      // Sort in descending order of y coordinates
      reverse(intersections_temp.begin(), intersections_temp.end());
    }

    // Converting to sorted vector of points
    for (y_sort point : intersections_temp) {
      intersections_sorted.push_back(Point(point.x, point.y));
    }
  }

  return intersections_sorted;

}
