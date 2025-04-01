#include "constants.hpp"
#include "inset_state.hpp"

// Function to automatically color topology based on contiguity graph
void InsetState::auto_color()
{
  // Colors are already provided
  if (colors_.size() == n_geo_divs())
    return;

  // Some colors are already provided, use light grey for the rest
  if (colors_.size() > 0) {
    for (GeoDiv &gd : geo_divs_) {
      if (!color_found(gd.id())) {
        insert_color(gd.id(), Color("#f2f2f2"));
      }
    }
  }

  std::vector<Color> palette;

  // Using default palette for now
  // TODO: Accept palette from user
  // From https://colorbrewer2.org/#type=qualitative&scheme=Pastel1&n=8
  palette.emplace_back("#fbb4ae");  // red
  palette.emplace_back("#b3cde3");  // blue
  palette.emplace_back("#ccebc5");  // green
  palette.emplace_back("#decbe4");  // purple
  palette.emplace_back("#fed9a6");  // orange
  palette.emplace_back("#ffffcc");  // yellow
  palette.emplace_back("#e5d8bd");  // brown
  palette.emplace_back("#fddaec");  // pink

  // Create continuity graph based on vertical and horizontal scanlines
  create_contiguity_graph();


  // Iterate until we are able to color the entire map
  for (const auto &gd : geo_divs_) {

    for (unsigned int i = 0; i < palette.size(); ++i) {
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
        // color_idx += color_idx_diff;
        break;
      }
    }
  }

  // Check if all GeoDivs are colored
  if (colors_.size() != n_geo_divs()) {
    std::cerr
      << "WARNING: Unable to color all GeoDivs without adjacency collisions"
      << std::endl;
    std::cerr
      << "Assigning #f2f2f2 (light grey) to uncolored GeoDivs"
      << std::endl;
      for (GeoDiv &gd : geo_divs_) {
        if (!color_found(gd.id())) {
          insert_color(gd.id(), Color("#f2f2f2"));
        }
      }
  }
}
