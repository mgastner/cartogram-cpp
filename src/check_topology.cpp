#include "cartogram_info.h"
#include "inset_state.h"
#include <CGAL/Boolean_set_operations_2.h>

// Returns error if there are holes not inside their respective polygons
void holes_inside_polygons(InsetState *inset_state)
{
  for (auto gd : inset_state->geo_divs()) {
    for (auto pwh : gd.polygons_with_holes()) {
      Polygon ext_ring = pwh.outer_boundary();
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        Polygon hole = *h;
        for (unsigned int i = 0; i < hole.size(); ++i) {

          // TODO: In the future, a better method would be to only check
          // whether one point in each hole is on the bounded side of the
          // exterior ring. Next, check whether the hole intersects the
          // polygon at any point. For this, the function
          // "do_intersect(Polygon, Polygon)" may help.
          if (ext_ring.bounded_side(hole[i]) == CGAL::ON_UNBOUNDED_SIDE) {
            CGAL::set_pretty_mode(std::cerr);
            std::cerr << "Hole detected outside polygon!" << std::endl;
            std::cerr << "Hole: " << hole << std::endl;
            std::cerr << "Polygon: " << ext_ring << std::endl;
            std::cerr << "GeoDiv: " << gd.id() << std::endl;
            _Exit(20);
          }
        }
      }
    }
  }
  return;
}


// Use machine epsilon (defined in constants.h) to get almost equal doubles.
// From https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
bool almost_equal(double a, double b) {
  return abs(a - b) <= dbl_epsilon * abs(a + b) * 2;
}

// Determine whether points are indistinguishable
bool points_almost_equal(Point a, Point b) {
  return (almost_equal(a[0], b[0]) && almost_equal(a[1], b[1]));
}

bool duplicates(std::vector<Point> v) {
  CGAL::set_pretty_mode(std::cerr);
  for (size_t i = 0; i < v.size() - 1; ++i) {
    if (points_almost_equal(v[i], v[i + 1])) {
      std::cerr << "i = " << i << std::endl;
      std::cerr << "Point: " << i << ", v[i]: " << v[i] << std::endl;
      std::cerr << "Point: "
                << i + 1
                << ", v[i + 1]: "
                << v[i + 1]
                << std::endl;
      return true;
    }
  }
  return false;
}
