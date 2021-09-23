#include "densification_points.h"

// This function takes in two lines a (defined by points a1 and a2) and b
// (defined by points b1 and b2) and returns the intersection between
// them, if any. If the two lines are parallel or are the same, return
// the point (-1, -1) that is always outside of the graticule grid cell.
XYPoint calc_intersection(XYPoint a1, XYPoint a2, XYPoint b1, XYPoint b2){
  
  // Get line equations
  double a = (a1.y - a2.y) / (a1.x - a2.x);
  double a_intercept = a1.y - a1.x * a;
  double b = (b1.y - b2.y) / (b1.x - b2.x);
  double b_intercept = b1.y - b1.x * b;

  XYPoint intersection;

  if (!isnan(a) && isnan(b)){

    // Case where only line b is vertical.
    intersection.x = b1.x;
    intersection.y = a * b1.x + a_intercept;

  } else if (isnan(a) && !isnan(b)){

    // Case where only line a is vertical.
    intersection.x = a1.x;
    intersection.y = b * a1.x + b_intercept;

  } else if ((!isnan(a) && !isnan(a)) || (a == b)){
    
    // No vertical line.
    intersection.x = (b_intercept - a_intercept) / (a - b);
    intersection.y = a * intersection.x + a_intercept;

  } else{
    
    // Set negative intersection coordinates if there is no solution or
    // infinitely many solutions.
    intersection.x = -1;
    intersection.y = -1;

  }

  return intersection;
}

// This function takes two points, "a" and "b", and returns all horizontal and
// vertical intersections with a graticule with graticule lines placed
// 1 unit apart. Further, it returns all points that collide with the diaganol
// of grid cells. Assumes grid cells to start at 0.5, 0.5
std::vector<Point> densification_points(Point a_, Point b_)
{
  
  // Vector for storing intersections before removing duplicates
  std::vector<XYPoint> temp_intersections;

  // Switch a and b if needed. [ab] and [ba] describe the same segment, but
  // if we flip the order of a and b, the resulting intersections will not
  // necessarily be the same due to floating point errors.
  XYPoint a;
  XYPoint b;
  if ((a_[0] > b_[0]) || ((a_[0] == b_[0]) && (a_[1] > b_[1]))){
    a.x = b_[0]; a.y = b_[1];
    b.x = a_[0]; b.y = a_[1];
  } else{
    a.x = a_[0]; a.y = a_[1];
    b.x = b_[0]; b.y = b_[1];
  }
  temp_intersections.push_back(a);
  temp_intersections.push_back(b);

  // Get graticule coordinates for a.
  XYPoint av0;
  av0.x = floor(a.x + 0.5) - 0.5;
  av0.y = floor(a.y + 0.5) - 0.5;

  // Get graticule coordinates for b.
  XYPoint bv0;
  bv0.x = floor(b.x + 0.5) - 0.5;
  bv0.y = floor(b.y + 0.5) - 0.5;

  // Get bottom-left (start_v) and top-right (end_v) corners of
  // the graticule cell rectangle (the smallest rectangular section
  // of the graticule grid cell containing both points)
  XYPoint start_v;
  XYPoint end_v;
  start_v.x = av0.x;
  end_v.x = bv0.x;
  if (a.y <= b.y){
    start_v.y = av0.y;
    end_v.y = bv0.y;
  } else {
    start_v.y = bv0.y;
    end_v.y = av0.y;
  }

  // Loop through each row, from bottom to top
  for (int i = int(start_v.y); i < int(end_v.y) + 1; ++i){
    
    // Loop through each column, from left to right
    for (int j = int(start_v.x); j < int(end_v.x) + 1; ++j){
      
      // Get points for the current graticule cell, in this order:
      // bottom-left, bottom-right, top-right, top-left
      XYPoint v0; v0.x = double(j) + 0.5; v0.y = double(i) + 0.5;
      XYPoint v1; v1.x = v0.x + 1.0; v1.y = v0.y;
      XYPoint v2; v2.x = v0.x + 1.0; v2.y = v0.y + 1.0;
      XYPoint v3; v3.x = v0.x; v3.y = v0.y + 1.0;

      std::vector<XYPoint> graticule_intersections;

      // Bottom intersection
      graticule_intersections.push_back(calc_intersection(a, b, v0, v1));
      // Left intersection
      graticule_intersections.push_back(calc_intersection(a, b, v0, v3));
      // Right intersection
      graticule_intersections.push_back(calc_intersection(a, b, v1, v2));
      // Top intersection 
      graticule_intersections.push_back(calc_intersection(a, b, v3, v2));
      // Diagonal intersections
      graticule_intersections.push_back(calc_intersection(a, b, v0, v2));
      graticule_intersections.push_back(calc_intersection(a, b, v3, v1));

      // Add segment intersections only
      for (XYPoint inter : graticule_intersections){
        if ((a.x <= inter.x && inter.x <= b.x) ||
            (b.x <= inter.x && inter.x <= a.x)){
          temp_intersections.push_back(inter);
        }
      }
    }
  }

  // Sort intersections
  std::sort(temp_intersections.begin(), temp_intersections.end());

  // Reverse if needed
  if ((a_[0] > b_[0]) || ((a_[0] == b_[0]) && (a_[1] > b_[1]))){
    std::reverse(temp_intersections.begin(), temp_intersections.end());
  }

  // Eliminate duplicates
  std::vector<Point> intersections;
  intersections.push_back(
    Point(temp_intersections[0].x, temp_intersections[0].y)
  );
  for (unsigned int i = 1; i < temp_intersections.size(); ++i){
    if ((temp_intersections[i - 1].x != temp_intersections[i].x) ||
        (temp_intersections[i - 1].y != temp_intersections[i].y)){
      intersections.push_back(
        Point(temp_intersections[i].x, temp_intersections[i].y)
      );
    }
  }
  return intersections;
}

// Round a double down to round_digits digits after decimal point.
double round_coordinate(double value){
  return round(value * round_digits) / round_digits;
}

// Truncate the coordinates of a point
Point round_point(Point p1) {
  double x1 = round_coordinate(p1[0]);
  double y1 = round_coordinate(p1[1]);
  Point rounded_point(x1, y1);
  return rounded_point;
}

// Very similar doubles
bool almost_equal(double a, double b) {
  return abs(a - b) <= (1/round_digits);
}

// Very similar points
bool point_almost_equal(Point a, Point b) {
  return (almost_equal(a[0], b[0]) && almost_equal(a[1], b[1]));
}

// // Returns ceiling up to nearest 0.5 value, e.g. 2.64 returns 3.5
// double half_ceil(double num) {

//   // Extracting decimal
//   double decimal = num - floor(num);

//   // Checking whether decimal is 0.5
//   if (abs(decimal - 0.5) <= 1e-10) {
//     decimal = 0.5;
//   }

//   // Checking whether to ceil
//   if (decimal > 0.5) {
//     return ceil(num) + 0.5;
//   }

//   // Else flooring
//   return floor(num) + 0.5;
// }

// // Returns floor up to nearest 0.5 value, e.g. 2.23 returns 1.5
// double half_floor(double num){
//   return floor(num + 0.5) - 0.5;
// }

// // Returns intersection of line with vertical grid line
// double calc_y_intersection(Point a, Point b, double x) {
//   return (a[1] * (b[0] - x) + b[1] * (x - a[0])) / (b[0] - a[0]);
// }

// // Returns intersection of line with horizontal grid line
// double calc_x_intersection(Point a, Point b, double y) {
//   return (a[0] * (b[1] - y) + b[0] * (y - a[1])) / (b[1] - a[1]);
// }

// // Returns intersection of line with vertical grid line
// double calc_y_intersection(XYPoint a, XYPoint b, double x) {
//   return (a.y * (b.x - x) + b.y * (x - a.x)) / (b.x - a.x);
// }

// // Returns intersection of line with horizontal grid line
// double calc_x_intersection(XYPoint a, XYPoint b, double y) {
//   return (a.x * (b.y - y) + b.x * (y - a.y)) / (b.y - a.y);
// }

// // This function takes two points, "a" and "b", and returns all horizontal and
// // vertical intersections with a graticule with graticule lines placed
// // 1 unit apart. Further, it returns all points that collide with the diaganol
// // of grid cells. Assumes grid cells to start at 0.5, 0.5
// std::vector<Point> densification_points(Point a, Point b)
// {

//   // Vector to store all intersections
//   std::vector<Point> intersections;
//   intersections.push_back(a);
//   intersections.push_back(b);

//   // Vertical intersections (intersections with x coordinate lines)
//   if (a[0] < b[0]) {

//     // X of "a" < x of "b", "x" is each vertical graticule line
//     for (double x = half_ceil(a[0]); x < b[0]; ++x) {
//       if (!almost_equal(x, a[0]) && !almost_equal(x, b[0])) {

//         // Get y coordinate, x coordinate = x
//         double y = calc_y_intersection(a, b, x);
//         Point temp(x, y);
//         if (!point_almost_equal(a, temp) && !point_almost_equal(b, temp)) {
//           intersections.push_back(round_point(temp));
//         }
//       }
//     }
//   } else if (a[0] > b[0]) {

//     // x of "a" > x of "b"
//     for (double x = half_floor(a[0]); x > b[0]; --x) {
//       if (!almost_equal(x, a[0]) && !almost_equal(x, b[0])) {

//         // Get y coordinate, x coordinate = x
//         double y = calc_y_intersection(a, b, x);
//         Point temp(x, y);
//         if (!point_almost_equal(a, temp) && !point_almost_equal(b, temp)) {
//           intersections.push_back(round_point(temp));
//         }
//       }
//     }
//   }

//   // Horizontal intersections (intersections with y coordinate lines)
//   if (a[1] < b[1]) {

//     // Y of "a" < y of "b", "y" is each horizontal graticule line
//     for (double y = half_ceil(a[1]); y < b[1]; ++y) {
//       if (!almost_equal(y, a[1]) && !almost_equal(y, b[1])) {
//         // Get x coordinate, y coordinate = y
//         double x = calc_x_intersection(a, b, y);
//         Point temp(x, y);

//         // Ensuring no corner points are pushed back as they would've already
//         // been added when calculating vertical intersections
//         if (x != half_floor(x) || a[0] == b[0]) {
//           if (!point_almost_equal(a, temp) && !point_almost_equal(b, temp)) {
//             intersections.push_back(round_point(temp));
//           }
//         }
//       }
//     }
//   } else if (a[1] > b[1]) {

//     // y of "a" > y of "b"
//     for (double y = half_floor(a[1]); y > b[1]; --y) {
//       if (!almost_equal(y, a[1]) && !almost_equal(y, b[1])) {

//         // Get x coordinate, y coordinate = y
//         double x = calc_x_intersection(a, b, y);
//         Point temp(x, y);
//         if (!point_almost_equal(a, temp) && !point_almost_equal(b, temp)) {
//           intersections.push_back(round_point(temp));
//         }
//       }
//     }
//   }

//   // Sorting intersections
//   std::sort(intersections.begin(), intersections.end());

//   // Clearing duplicates
//   intersections.erase(unique(intersections.begin(), intersections.end()),
//   intersections.end());

//   // Storing current size
//   double size = intersections.size();

//   // TODO: Fix bug for diaganols when 0 intersections

//   // Inserting intersection with diagonals

//   // Adding intersections with bottom_right-top_left diaganol
//   for (unsigned int i = 0; i < size - 1; ++i) {

//     // std::cout << "Inside for loop! Calculating diaganol intersections...\n";;

//     // Bottom-right corner
//     Point bottom_right(half_floor(intersections[i][0]) + 1,
//                        half_floor(intersections[i][1]));

//     // Top-left corner
//     Point top_left(bottom_right[0] - 1,
//                    bottom_right[1] + 1);

//     // Finding intersection
//     Segment seg_diaganol(top_left, bottom_right);
//     Segment seg_intersec(a, b);
//     auto result = CGAL::intersection(seg_diaganol, seg_intersec);
//     if (result) {
//       if (const Point *s = boost::get<Point>(&*result)) {
//         if (!point_almost_equal(a, (*s)) && !point_almost_equal(b, (*s))) {
//           intersections.push_back(round_point((*s)));
//         }
//       }
//     }
//   }

//   // Adding intersections with bottom_left-top_right diaganol
//   for (unsigned int i = 0; i < size - 1; ++i) {

//     // Bottom-left corner
//     Point bottom_left(half_floor(intersections[i][0]),
//                       half_floor(intersections[i][1]));

//     // Top-right corner
//     Point top_right(bottom_left[0] + 1,
//                     bottom_left[1] + 1);

//     // Finding intersection
//     Segment seg_diaganol(top_right, bottom_left);
//     Segment seg_intersec(a, b);

//     auto result = CGAL::intersection(seg_diaganol, seg_intersec);
//     if (result) {
//       if (const Point *s = boost::get<Point>(&*result)) {

//         // Checking intersection not through middle of graticule lines,
//         // Because if middle of graticule lines, then duplicates will be added as
//         // intersection already added by previous for loop.
//         if (!almost_equal((*s)[0], floor((*s)[0]))) {
//           if (!point_almost_equal(a, (*s)) && !point_almost_equal(b, (*s))) {
//             intersections.push_back(round_point((*s)));
//           }
//         }
//       }
//     }
//   }

//   // TODO: REMOVE DUPLICATES UNTIL NO MORE

//   // Sorting intersections
//   std::sort(intersections.begin(), intersections.end());

//   // Removing duplicates
//   intersections.erase(unique(intersections.begin(), intersections.end()),
//   intersections.end());

//   if (a[0] > b[0] || (a[1] > b[1] && almost_equal(a[0], b[0]))) {

//     // Sort in descending order of x coordinates
//     std::reverse(intersections.begin(), intersections.end());
//   }
  
//   return intersections;
// }
