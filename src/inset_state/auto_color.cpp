#include "../constants.h"
#include "../cartogram_info.h"
#include "../inset_state.h"
#include "../colors.h"

// Function to automatically color topology based on contiguity graph
void InsetState::auto_color()
{

  std::vector<Color> palette;

  // Using default Palette for now
  // From https://colorbrewer2.org/
  palette.push_back(Color(27, 158, 119)); // turquoise
  palette.push_back(Color(117, 112, 179)); // purple
  palette.push_back(Color(231, 41, 138)); // pinkish-magenta
  palette.push_back(Color(102, 166, 30)); // green
  palette.push_back(Color(230, 171, 2)); // yellow
  palette.push_back(Color(217, 95, 2)); // redish-orange
  palette.push_back(Color(102, 102, 102)); // dark grey
  palette.push_back(Color(166, 118, 29)); // mustard-brown

  // Find resolution
  unsigned int res = default_res;

  // Creating full continuity graph based on vertical and horizontal scanlines
  create_contiguity_graph(res);

  // Count to maximize colors used. This changes the starting color that
  // the algorithm choses for each GeoDiv.
  unsigned int count = 0;

  // If coloring was unsucesful, we restrict ourselves to a smaller count.
  // We begin by using as many colors as possible.
  int max_i = palette.size();

  // Iterating until we are able to color the entire map
  while (colors_size() < n_geo_divs() &&
         max_i >= 0) {

    // Iterating through GeoDivs
    for (auto gd : geo_divs()) {

      // Iterating through all possible colors
      for (size_t i = (count % max_i); i < palette.size(); ++i) {
        Color c = palette[i];
        bool shared_color = false;

        // Iterating through adjacent GeoDivs to check whether shared color
        for (std::string gd_id : gd.adjacent_geodivs()) {

          // Checking to see whether color found
          if (color_found(gd_id)) {
            if (colors_at(gd_id) == c) {
              shared_color = true;
            }
          }
        }

        // If the color is not shared with any other adjacent GeoDiv,
        // we assign it the color.
        if (!shared_color) {
          colors_insert(gd.id(), c);
          i = palette.size();
        }
      }
      count++;
    }
    max_i--;
  }
  return;
}
