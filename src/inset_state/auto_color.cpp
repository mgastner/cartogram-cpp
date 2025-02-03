#include "constants.hpp"
#include "inset_state.hpp"

// Function to automatically color topology based on contiguity graph
void InsetState::auto_color()
{
  // Colors are already provided
  if (!colors_.empty()) return;

  std::vector<Color> palette;

  // Using default palette for now
  // TODO: Accept palette from user
  // From https://colorbrewer2.org/#type=qualitative&scheme=Set3&n=8
  palette.emplace_back("#8dd3c7");  // green
  palette.emplace_back("#ffffb3");  // yellow
  palette.emplace_back("#bebada");  // purple
  palette.emplace_back("#fb8072");  // red
  palette.emplace_back("#80b1d3");  // blue
  palette.emplace_back("#fdb462");  // orange
  palette.emplace_back("#b3de69");  // green
  palette.emplace_back("#fccde5");  // pink

  // Find resolution
  const unsigned int resolution = intersections_resolution;

  // Create continuity graph based on vertical and horizontal scanlines
  create_contiguity_graph(resolution);

  // Count colors used
  unsigned int count = 0;

  // If coloring was unsuccessful, we restrict ourselves to a smaller count.
  // We begin by using as many colors as possible.
  int max_i = palette.size();

  // Iterate until we are able to color the entire map
  while (colors_.size() < n_geo_divs() && max_i >= 0) {
    for (const auto &gd : geo_divs_) {

      // Iterate over all possible colors
      for (unsigned int i = (count % max_i); i < palette.size(); ++i) {
        const Color c = palette[i];
        bool shared_color = false;

        // Check whether adjacent GeoDivs have the same color
        for (const auto &gd_id : gd.adjacent_geodivs()) {
          if (color_found(gd_id)) {
            if (color_at(gd_id) == c) {
              shared_color = true;
            }
          }
        }

        // Assign color if it is not shared with any adjacent GeoDiv
        if (!shared_color) {
          insert_color(gd.id(), c);
          i = palette.size();
        }
      }
      count++;
    }
    max_i--;
  }
}
