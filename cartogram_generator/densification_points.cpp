#include "map_state.h"
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

  //return floor(num + 0.5) + 0.5;

}

double simple_half_ceil(double num){
  return floor(num + 0.5) + 0.5;
}

// Returns floor up to nearest 0.5 value, e.g. 2.23 returns 1.5
double half_floor(double num) {

  // Extracting decimal
  double decimal = num - floor(num);

  // Checking whether decimal is 0.5
  if (abs(decimal - 0.5) < 1e-10) {
    decimal = 0.5;
  }

  // Checking whether to add 0.5
  if (decimal >= 0.5) {
    return floor(num) + 0.5;
  }

  // Else subtracting 0.5
  return floor(num) - 0.5;

  //return floor(num + 0.5) - 0.5;
}

double simple_half_floor(double num){
  return floor(num + 0.5) - 0.5;
}

bool almost_equal(double a, double b) {

  // Very similar doubles
  return abs(a - b) <= 1e-10;
}

bool point_almost_equal(Point a, Point b) {

  // Very similar points
  return (almost_equal(a[0], b[0]) && almost_equal(a[1], b[1]));
}

double round_down(double value) {
    return floor(value * 1e10) / 1e10;
}

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

  // CGAL::set_pretty_mode(std::cout);
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

  // std::cout << "Vertical intersections added!" << '\n';
  // std::cout << "Current state: " << '\n';
  // for (size_t i = 0; i < intersections.size(); i++) {
  //   std::cout << i + 1 << ": ";
  //   std::cout << intersections[i] << '\n';
  // }

  // Horizontal intersections (intersections with y coordinate lines)
  if (a[1] < b[1]) {

    // Y of "a" < y of "b", "y" is each horizontal graticule line
    for (double y = half_ceil(a[1]); y < b[1]; ++y) {
      if (!almost_equal(y, a[1]) && !almost_equal(y, b[1])) {
        // Get x coordinate, y coordinate = y
        double x = calc_x_intersection(a, b, y);
        Point temp(x, y);
        // double difference = 3.5 - x;
        // printf("x with printf: %.10e\n", x);
        // printf("difference with printf: %.10e\n", difference);

        // Ensuring no corner points are pushed back
        // As they would've already been added when calculating
        // Vertical intersections
        // if ((!almost_equal(x, half_floor(x)) && !almost_equal(y, half_floor(y))) || a[0] == b[0]) {
        if (x != half_floor(x) || a[0] == b[0]) {

          // std::cout << "x != half_floor(x) ?" << '\n';
          // x != half_floor(x) ? std::cout << "true" : std::cout << "false";
          // std::cout << '\n';
          // std::cout << "x: " << std::setprecision(10) << x << '\n';
          // std::cout << "half_floor(x): " << half_floor(x) << '\n';
          // std::cout << "a[0] == b[0] ?" << '\n';
          // a[0] == b[0] ? std::cout << "true" : std::cout << "false";
          // std::cout << '\n';
          if (!point_almost_equal(a, temp) && !point_almost_equal(b, temp)) {
            intersections.push_back(round_point(temp));
          }
        }
      }
    }
  } else if (a[1] > b[1]) {

    // y of "a" > y of "b"
    // for (double y = half_floor(a[1]); y > b[1]; --y) {
    //   if (!almost_equal(y, a[1]) && !almost_equal(y, b[1])) {

    //     // Get x coordinate, y coordinate = y
    //     double x = calc_x_intersection(a, b, y);
    //     Point temp(x, y);
    //     // Ensuring no corner points are pushed back
    //     if ((!almost_equal(x, half_floor(x)) && !almost_equal(y, half_floor(y))) || a[0] == b[0]) {
    //       if (!point_almost_equal(a, temp) && !point_almost_equal(b, temp)) {
    //         intersections.push_back(round_point(temp));
    //       }
    //     }
    //   }
    // }
    for (double y = half_floor(a[1]); y > b[1]; --y) {
      if (!almost_equal(y, a[1]) && !almost_equal(y, b[1])) {

        // Get x coordinate, y coordinate = y
        double x = calc_x_intersection(a, b, y);
        Point temp(x, y);
        // Ensuring no corner points are pushed back
        // if ((!almost_equal(x, half_floor(x)) && !almost_equal(y, half_floor(y))) || a[0] == b[0]) {
        //   if (!point_almost_equal(a, temp) && !point_almost_equal(b, temp)) {
        //     intersections.push_back(round_point(temp));
        //   }
        // }
        if (!point_almost_equal(a, temp) && !point_almost_equal(b, temp)) {
          intersections.push_back(round_point(temp));
        }
      }
    }
  }

  // std::cout << "Horizontal intersections added!" << '\n';
  // std::cout << "Current state: " << '\n';
  // for (size_t i = 0; i < intersections.size(); i++) {
  //   std::cout << i + 1 << ": ";
  //   std::cout << intersections[i] << '\n';
  // }

  // std::cout << "Graticule intersections calculated; Now sorting..." << '\n';

  // Sorting intersections
  std::sort(intersections.begin(), intersections.end());

  // std::cout << "State after sorting:" << '\n';
  // for (size_t i = 0; i < intersections.size(); i++) {
  //   std::cout << i + 1 << ": ";
  //   std::cout << intersections[i] << '\n';
  // }

  // Clearing duplicates
  intersections.erase(unique(intersections.begin(), intersections.end()),
  intersections.end());

  // std::cout << "Element 2: " << intersections[1] << '\n';
  // std::cout << "Element 3: " << intersections[2] << '\n';
  // std::cout << "Element 3 == Element 2 ? ";
  // intersections[1] == intersections[2] ? std::cout << "yes" : std::cout << "no";
  // std::cout << '\n';

  // std::cout << "State after clearing:" << '\n';
  // for (size_t i = 0; i < intersections.size(); i++) {
  //   std::cout << i + 1 << ": ";
  //   std::cout << intersections[i] << '\n';
  // }

  // // actual size
  // std::set<Point> counter(intersections.begin(), intersections.end());
  // intersections.assign(intersections.begin(), intersections.end());
  //
  // // for (int i = 0; i < intersections.size(); i++)
  // //       counter.insert(intersections[i]);
  // // double actual_size = counter.size();
  //
  // double actual_size = counter.size();
  //
  // std::cout << "Actual size: " << actual_size << '\n';
  // std::cout << "Intersections size: " << intersections.size() << '\n';
  //
  // while (actual_size < intersections.size()) {
  //   std::cout << "Does this even run?" << '\n';
  //   intersections.erase(unique(intersections.begin(), intersections.end()),
  //   intersections.end());
  // }
  //
  //
  // std::cout << "Duplicates cleared!" << '\n';
  // std::cout << "Current state: " << '\n';
  // for (size_t i = 0; i < intersections.size(); i++) {
  //   std::cout << i + 1 << ": ";
  //   std::cout << intersections[i] << '\n';
  // }

  // Storing current size
  double size = intersections.size();

  // TODO: Fix bug for diaganols when 0 intersections
  // TODO: Round to 10 decimal places

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
    Segment seg_intersec(a, b);

    auto result = intersection(seg_diaganol, seg_intersec);
    if (result) {
      if (const Point *s = boost::get<Point>(&*result)) {
        if (!point_almost_equal(a, (*s)) && !point_almost_equal(b, (*s))) {
          intersections.push_back(round_point((*s)));
        }
      }
    }
  }

  // std::cout << "bottom_right-top_left intersections added!" << '\n';
  // std::cout << "Current state: " << '\n';
  // for (size_t i = 0; i < intersections.size(); i++) {
  //   std::cout << i + 1 << ": ";
  //   std::cout << intersections[i] << '\n';
  // }

  // Adding intersections with bottom_left-top_right diaganol
  for (unsigned int i = 0; i < size - 1; ++i) {

    // Bottom-left corner
    Point bottom_left(half_floor(intersections[i][0]),
                      half_floor(intersections[i][1]));

    // Top-right corner
    Point top_right(half_ceil(intersections[i + 1][0]),
                    half_ceil(intersections[i + 1][1]));

    if (a[0] >= 197.485 && a[0] <= 197.486 && a[1] >= 257.620 && a[1] <= 257.621 &&
        b[0] >= 197.667 && b[0] <= 197.668 && b[1] >= 257.403 && b[1] <= 257.404) {
      
      std::cout << "\nCurrent diagonal: \n"
                << "(" << bottom_left[0] << ", " << bottom_left[1] << ")\n"
                << "(" << top_right[0] << ", " << top_right[1] << ")\n";

      std::cout << "Rounded from: \n"
                << "(" << intersections[i][0] << ", " << intersections[i][1] << ")\n"
                << "(" << intersections[i + 1][0] << ", " << intersections[i + 1][1] << ")\n";
    }

    // Finding intersection
    Segment seg_diaganol(top_right, bottom_left);
    Segment seg_intersec(a, b);

    auto result = intersection(seg_diaganol, seg_intersec);
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

  // std::cout << "bottom_left-top_right intersections added!" << '\n';
  // std::cout << "Current state: " << '\n';
  // for (size_t i = 0; i < intersections.size(); i++) {
  //   std::cout << i + 1 << ": ";
  //   std::cout << intersections[i] << '\n';
  // }

  // if (size == 0) {
  //
  //   std::cout << "if size 0, intersections added!" << '\n';
  //   // Adding intersection with singular grid cell diaganols
  //
  //   // Bottom-right corner
  //   Point bottom_right(half_floor(a[0]) + 1,
  //                      half_floor(a[1]));
  //
  //   // Top-left corner
  //   Point top_left(half_ceil(a[0]) - 1,
  //                  half_ceil(a[1]));
  //
  //   // Bottom-left corner
  //   Point bottom_left(half_floor(a[0]),
  //                     half_floor(a[1]));
  //
  //   // Top-right corner
  //   Point top_right(half_ceil(a[0]),
  //                   half_ceil(a[1]));
  //
  //   // Finding intersections
  //   Segment seg_tlbr(top_left, bottom_right);
  //   Segment seg_trbl(top_right, bottom_left);
  //   Segment seg_intersec(a, b);
  //
  //   auto result_1 = intersection(seg_tlbr, seg_intersec);
  //   if (result_1) {
  //     if (const Point *s = boost::get<Point>(&*result_1)) {
  //       intersections.push_back(round_point(*s));
  //     }
  //   }
  //
  //   auto result_2 = intersection(seg_trbl, seg_intersec);
  //   if (result_2) {
  //     if (const Point *s = boost::get<Point>(&*result_2)) {
  //       if (!((*s)[0] == floor((*s)[0]) || (*s)[1] == floor((*s)[1]))) {
  //         intersections.push_back(round_point(*s));
  //       }
  //     }
  //   }

  // std::cout << "Current state: " << '\n';
  // for (size_t i = 0; i < intersections.size(); i++) {
  //   std::cout << i + 1 << ": ";
  //   std::cout << intersections[i] << '\n';
  // }
  //
  // } else {
  //
  //   std::cout << "size > 0, final intersection added!" << '\n';
  //
  //   // Adding final intersection using last intersection point
  //
  //   // Bottom-right corner
  //   Point bottom_right(half_floor(a[0]) + 1,
  //                      half_floor(a[1]));
  //
  //   // Top-left corner
  //   Point top_left(half_ceil(a[0]) - 1,
  //                  half_ceil(a[1]));
  //
  //   // Bottom-left corner
  //   Point bottom_left(half_floor(a[0]),
  //                     half_floor(a[1]));
  //
  //   // Top-right corner
  //   Point top_right(half_ceil(a[0]),
  //                   half_ceil(a[1]));
  //
  //   // Finding intersections
  //   Segment seg_tlbr(top_left, bottom_right);
  //   Segment seg_trbl(top_right, bottom_left);
  //
  //   Point c;
  //   if (a[0] > b[0] || (a[1] > b[1] && a[0] == b[0])) {
  //     c = a;
  //   } else {
  //     c = b;
  //   }
  //
  //   Segment seg_intersec(intersections[size - 1], c);
  //
  //   auto result_1 = intersection(seg_tlbr, seg_intersec);
  //   if (result_1) {
  //     if (const Point *s = boost::get<Point>(&*result_1)) {
  //       intersections.push_back(round_point(*s));
  //     }
  //   }
  //
  //   auto result_2 = intersection(seg_trbl, seg_intersec);
  //   if (result_2) {
  //     if (const Point *s = boost::get<Point>(&*result_2)) {
  //       if (!((*s)[0] == floor((*s)[0]) || (*s)[1] == floor((*s)[1]))) {
  //         intersections.push_back(round_point(*s));
  //       }
  //     }
  //   }
  // }

  // TODO: REMOVE DUPLICATES UNTIL NO MORE

  // std::cout << "Diaganol intersections calculated! Now sorting..." << '\n';

  // Sorting intersections
  std::sort(intersections.begin(), intersections.end());

  // Removing duplicates
  intersections.erase(unique(intersections.begin(), intersections.end()),
  intersections.end());

  if (a[0] > b[0] || (a[1] > b[1] && almost_equal(a[0], b[0]))) {

    // Sort in descending order of x coordinates
    std::reverse(intersections.begin(), intersections.end());
  }

  // std::cout << "Sorted! Now returning..." << '\n';

  return intersections;
}
