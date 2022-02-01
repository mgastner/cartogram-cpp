#ifndef XYPOINT_H_
#define XYPOINT_H_

// Struct to store the X and Y coordinates of a 2D point as doubles.
// Calculations with double coordinates produce more predictable
// floating-point truncation than CGAL's Point data structure.
struct XYPoint {
  double x;
  double y;

  // Flip x and y coordinates
  void flip() {
    double temp = x;
    x = y;
    y = temp;
  }

  // Overload "<" operator. Idea from
  // https://stackoverflow.com/questions/4892680/sorting-a-vector-of-structs
  bool operator < (const XYPoint &rhs) const
  {
    if (x != rhs.x) {
      return (x < rhs.x);
    }
    return y < rhs.y;
  }

  // Overload "==" operator
  bool operator == (const XYPoint &rhs) const
  {
    return (x == rhs.x && y == rhs.y);
  }

  // Constructor, setting points to 0
  XYPoint()
  {
    x = 0;
    y = 0;
  }

  // Constructor with two values given
  XYPoint(double xg, double yg) {
    x = xg;
    y = yg;
  }
};

#endif
