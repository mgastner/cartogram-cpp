#include "../inset_state.h"
#include "../constants.h"

void InsetState::create_contiguity_graph()
{
  // Assign integers to GeoDivs
  std::map<std::string, unsigned int> int_at_gd;
  size_t i = 0;
  for (const auto &gd : geo_divs()) {
    int_at_gd[gd.id()] = i;
    i++;
  }

  // For all points in the inset, determine which GeoDivs have this point in
  // their boundaries
  std::map<Point, std::vector<unsigned int> > gds_at_point;
  for (const auto &gd : geo_divs()) {
    const auto gd_as_int = int_at_gd.at(gd.id());
    for (const auto &pwh : gd.polygons_with_holes()) {
      for (const auto &p : pwh.outer_boundary()) {
        gds_at_point[p].push_back(gd_as_int);
      }
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        for (const auto &p : (*h)) {
          gds_at_point[p].push_back(gd_as_int);
        }
      }
    }
  }

  // Iterate over all points
  for (const auto &point_info : gds_at_point) {
    const auto &geo_divs_with_point = point_info.second;
    for (unsigned int i = 0; i < geo_divs_with_point.size(); ++i) {
      for (unsigned int j = i + 1; j < geo_divs_with_point.size(); ++j) {
        const auto gd_i = geo_divs_with_point[i];
        const auto gd_j = geo_divs_with_point[j];

        // Ensure that GeoDivs are different
        if (gd_i != gd_j) {
          geo_divs_[gd_i].adjacent_to(geo_divs_[gd_j].id());
          geo_divs_[gd_j].adjacent_to(geo_divs_[gd_i].id());
        }
      }
    }
  }
}
