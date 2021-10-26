#include "densification_points.h"

// Use machine epsilon (defined in constants.h) to get almost equal doubles.
// From https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
bool almost_equal(double a, double b) {
  return fabs(a - b) <= dbl_epsilon * fabs(a + b) * 2;
}

// Very similar points
bool point_almost_equal(Point a, Point b) {
  return (almost_equal(a[0], b[0]) && almost_equal(a[1], b[1]));
}


// This function takes in two lines a (defined by points a1 and a2) and b
// (defined by points b1 and b2) and returns the intersection between
// them, if any. If the two lines are parallel or are the same, return
// the point (-1, -1) that is always outside of the graticule grid cell.
XYPoint calc_intersection(XYPoint a1, XYPoint a2, XYPoint b1, XYPoint b2){
  
  // Check if any segment is undefined, i.e. defined by identical points
  if (a1 == a2 || b1 == b2) {
    std::cerr << "ERROR: End points of line segment are identical" << std::endl;
    _Exit(EXIT_FAILURE);
  }
  
  // Get line equations
  double a = (a1.y - a2.y) / (a1.x - a2.x);
  double a_intercept = a1.y - (a1.x * a);
  double b = (b1.y - b2.y) / (b1.x - b2.x);
  double b_intercept = b1.y - (b1.x * b);

  XYPoint intersection;

  if (isfinite(a) && isfinite(b) && a != b){

    // No vertical line.
    intersection.x = (b_intercept - a_intercept) / (a - b);
    intersection.y = a * intersection.x + a_intercept;

  } else if (isfinite(a) && isinf(b)){

    // Case where only line b is vertical.
    intersection.x = b1.x;
    intersection.y = a * b1.x + a_intercept;

  } else if (isfinite(b) && isinf(a)){

    // Case where only line a is vertical.
    intersection.x = a1.x;
    intersection.y = b * a1.x + b_intercept;

  } else {

    // Set negative intersection coordinates if there is no solution or
    // infinitely many solutions.
    intersection.x = -1;
    intersection.y = -1;

  }

  return intersection;
}

// This function takes two points, "a" and "b", and returns all horizontal and
// vertical intersections with a graticule with graticule lines placed
// 1 unit apart. It also returns all intersections with the diagonals of these
// graticule cells. Assumes graticule cells to start at 0.5, 0.5.
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

  // Get bottom-left point of graticule cell containing a.
  XYPoint av0;
  av0.x = floor(a.x + 0.5) - 0.5;
  av0.y = floor(a.y + 0.5) - 0.5;

  // Get bottom-left point of graticule cell containing b.
  XYPoint bv0;
  bv0.x = floor(b.x + 0.5) - 0.5;
  bv0.y = floor(b.y + 0.5) - 0.5;

  // Get bottom-left (start_v) and top-right (end_v) graticule cells of
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

      // Add segment intersections only. Usually, it is enough to check that the x
      // coordinate of the intersection is within bounds, but in some edge cases,
      // it is possible that the x-coordinate is without bounds but the y coordinate
      // is outside, i.e. when (ab) is vertical.
      for (XYPoint inter : graticule_intersections){
        if (((a.x <= inter.x && inter.x <= b.x) || (b.x <= inter.x && inter.x <= a.x)) &&
            ((a.y <= inter.y && inter.y <= b.y) || (b.y <= inter.y && inter.y <= a.y))){
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
    if (temp_intersections[i - 1] != temp_intersections[i]){
      intersections.push_back(
        Point(temp_intersections[i].x, temp_intersections[i].y)
      );
    }
  }
  return intersections;
}
