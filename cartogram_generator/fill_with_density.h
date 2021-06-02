#ifndef FILL_WITH_DENSITY_H_
#define FILL_WITH_DENSITY_H_

#include <string.h>

// Struct to store intersection data
struct intersection {
  double x;  // x-coordinate of intersection
  double target_density;  // GeoDiv's target_density
  std::string geo_div_id; // GeoDIv's ID
  bool direction;  // Does intersection enter (true) or exit (false)?

  // Overload "<" operator for this data type. Idea from
  // https://stackoverflow.com/questions/4892680/sorting-a-vector-of-structs
  bool operator < (const intersection &rhs) const
  {
    return (x < rhs.x || (x == rhs.x && direction < rhs.direction));
  }
};

// Contribution to a cell's density that comes from one ray and one geo_div
struct density { // density_contribution
  double weight; // Store weight (area_err * how much of ray is inside geo_div)
  double target_density; // Store target density

  density ()
  {
    weight = 0;
    target_density = 0;
  }

  density (double weight_, double target_density_)
  {
    weight = weight_;
    target_density = target_density_;
  }

};

// Vector of all density contributions to a particular cell
typedef std::vector<density> cell; // cell_dens_info

void fill_with_density(MapState*);

#endif
