#include "cartogram_info.h"
#include "colors.h"
#include "constants.h"
#include "inset_state.h"

// Function to automatically color topology based on contiguity graph
void InsetState::auto_color()
{
  std::vector<Color> palette;

  // Using default palette for now
  // TODO: Accept palette from user
  // From https://colorbrewer2.org/#type=qualitative&scheme=Set3&n=8
  // palette.push_back(Color("#8dd3c7"));  // green
  // palette.push_back(Color("#ffffb3"));  // yellow
  // palette.push_back(Color("#bebada"));  // purple
  // palette.push_back(Color("#fb8072"));  // red
  // palette.push_back(Color("#80b1d3"));  // blue
  // palette.push_back(Color("#fdb462"));  // orange
  // palette.push_back(Color("#b3de69"));  // green
  // palette.push_back(Color("#fccde5"));  // pink

  // Second color bar
  // palette.push_back(Color("#a6cee3"));  // light blue
  // palette.push_back(Color("#b2df8a"));  // light green
  // palette.push_back(Color("#fb9a99"));  // light red
  // palette.push_back(Color("#fdbf6f"));  // light orange
  // palette.push_back(Color("#cab2d6"));  // light purple
  // palette.push_back(Color("#33a02c"));  // dark green
  // palette.push_back(Color("#1f78b4"));  // dark blue
  // palette.push_back(Color("#e31a1c"));  // dark red
  // palette.push_back(Color("#ff7f00"));  // dark orange
  // palette.push_back(Color("#6a3d9a"));  // dark purple
  // palette.push_back(Color("#ffff99"));  // yellow
  // palette.push_back(Color("#b15928"));  // brown

  // From https://flatuicolors.com/palette/defo
  palette.push_back(Color("#e74c3c"));  // red
  palette.push_back(Color("#f1c40f"));  // yellow
  palette.push_back(Color("#34495e"));  // dark blue
  palette.push_back(Color("#e67e22"));  // orange
  palette.push_back(Color("#9b59b6"));  // amethyst
  palette.push_back(Color("#2ecc71"));  // green
  palette.push_back(Color("#3498db"));  // light blue
  palette.push_back(Color("#1abc9c"));  // turqiose

  // Find resolution
  const unsigned int resolution = default_resolution;

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
  return;
}
