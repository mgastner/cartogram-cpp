#ifndef FILL_WITH_DENSITY_H_
#define FILL_WITH_DENSITY_H_

// Struct to store intersection data
struct intersection {
  double x;  // x-coordinate of intersection
  double target_density;  // GeoDiv's target_density
  bool direction;  // Does intersection enter or exit?

  // Overload "<" operator for this data type. Idea from
  // https://stackoverflow.com/questions/4892680/sorting-a-vector-of-structs
  bool operator < (const intersection &rhs) const
  {
    return (x < rhs.x || (x == rhs.x && direction < rhs.direction));
  }
};

void fill_with_density(MapState*);

#endif
