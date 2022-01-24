#include "../inset_state.h"
#include "../constants.h"

void InsetState::create_contiguity_graph()
{

  // Finding positions of GeoDivs
  std::map<std::string, unsigned int> int_at_gd;
  size_t i = 0;
  for (auto gd : geo_divs()) {
    int_at_gd[gd.id()] = i;
    i++;
  }

  // Mapping points to geodivs (which GeoDivs have this point)
  std::map<Point, std::vector<unsigned int>> gds_at_point;

  for (auto gd : geo_divs()) {
    unsigned int gd_as_int = int_at_gd.at(gd.id());
    for (const auto &pwh : gd.polygons_with_holes()) {
      for (Point p : pwh.outer_boundary()) {
        gds_at_point[p].push_back(gd_as_int);
      }
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        for (Point p : (*h)) {
          gds_at_point[p].push_back(gd_as_int);
        }
      }
    }
  }

  // Iterating through all points
  for (auto [Point, geo_divs_with_point] : gds_at_point ) {
    for (unsigned int i = 0; i < geo_divs_with_point.size(); ++i) {
      for (unsigned int j = i + 1; j < geo_divs_with_point.size(); ++j) {
        unsigned int gd_i = geo_divs_with_point[i];
        unsigned int gd_j = geo_divs_with_point[j];

        // Ensuring that the geodivs are different
        if (gd_i != gd_j) {
          geo_divs_[gd_i].adjacent_to(geo_divs_[gd_j].id());
          geo_divs_[gd_j].adjacent_to(geo_divs_[gd_i].id());
        }
      }
    }
  }
}
