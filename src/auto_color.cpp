#include "cartogram_info.h"
#include "inset_state.h"
#include "auto_color.h"
#include "colors.h"

void create_adjacency_graph(InsetState* inset_state,
                            char graph,
                            unsigned int res)
{

  // Getting the chosen graph.
  std::vector<std::vector<intersection> > adj_graph;
  unsigned int max_k = 0;
  if (graph == 'h') {
    adj_graph = inset_state->horizontal_adj();
    max_k = inset_state->ly();
  } else if (graph == 'v') {
    adj_graph = inset_state->vertical_adj();
    max_k = inset_state->lx();
  }

  // Iterating through horizontal adjacency graph
  for (unsigned int k = 0; k < max_k; ++k) {

    // Cycle through each of the "res" number of rays in one cell
    for (double ray = k + (1.0/res)/2;
         ray < k + 1;
         ray += (1.0/res)) {

      // The intersections for one ray
      std::vector<intersection> intersections =
        adj_graph[static_cast<int>(round(((ray - (1.0/res)/2.0) * res)))];

      // Sort vector in ascending order of intersection
      sort(intersections.begin(), intersections.end());

      int size = intersections.size() - 1;

      // Fill GeoDivs by iterating through intersections
      for (int l = 1; l < size; l += 2) {
        double x_1 = intersections[l].coord;
        double x_2 = intersections[l + 1].coord;
        std::string gd_1 = intersections[l].geo_div_id;
        std::string gd_2 = intersections[l + 1].geo_div_id;

        if (gd_1 != gd_2 && x_1 == x_2) {
          for (auto &gd : *inset_state->ref_to_geo_divs()) {
            if (gd.id() == gd_1) {
              gd.adjacent_to(gd_2);
            } else if (gd.id() == gd_2) {
              gd.adjacent_to(gd_1);
            }
          }
        }

      }
    }
  }
}

// Function to automatically color topology based on adjacency graph
void auto_color(InsetState* inset_state)
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


  // Getting the horizontal adjacency graph.
  std::vector<std::vector<intersection> > horizontal_adj =
    inset_state->horizontal_adj();

  // Find resolution
  unsigned int res = horizontal_adj.size() / inset_state->ly();

  // Creating vertical adjacency graph
  inset_state->create_vertical_adjacency_graph(res);

  // Creating full adjacency graph based on vertical and horizontal graphs
  create_adjacency_graph(inset_state, 'h', res);
  create_adjacency_graph(inset_state, 'v', res);

  // Count to maximize colors used. This changes the starting color that
  // the algorithm choses for each GeoDiv.
  unsigned int count = 0;

  // If coloring was unsucesful, we restrict ourselves to a smaller count.
  // We begin by using as many colors as possible.
  int max_i = palette.size();

  // Iterating until we are able to color the entire map
  while (inset_state->colors_size() < inset_state->n_geo_divs() &&
         max_i >= 0) {

    // Iterating through GeoDivs
    for (auto &gd : *inset_state->ref_to_geo_divs()) {

      // Iterating through all possible colors
      for (size_t i = (count % max_i); i < palette.size(); ++i) {
        Color c = palette[i];
        bool shared_color = false;

        // Iterating through adjacent GeoDivs to check whether shared color
        for (std::string gd_id : gd.adjacent_geodivs()) {

          // Checking to see whether color found
          if (inset_state->color_found(gd_id)) {
            if (inset_state->colors_at(gd_id) == c) {
              shared_color = true;
            }
          }
        }

        // If the color is not shared with any other adjacent GeoDiv,
        // we assign it the color.
        if (!shared_color) {
          inset_state->colors_insert(gd.id(), c);
          i = palette.size();
        }
      }
      count++;
    }
    max_i--;
  }
  return;
}
