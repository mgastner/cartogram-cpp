#ifndef XY_POINT_H_
#define XY_POINT_H_

// Struct to store the X and Y coordinates of a 2D point as doubles.
// Calculations with double coordinates produce more predictable
// floating-point truncation than CGAL's Point data structure.
struct XYPoint {
  double x;
  double y;

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
};

#endif
