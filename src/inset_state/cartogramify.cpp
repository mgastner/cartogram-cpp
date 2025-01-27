#include "inset_state.hpp"

void InsetState::preprocess()
{
    // Remove tiny polygons below threshold
    if (args_.remove_tiny_polygons) {
        remove_tiny_polygons(args_.min_polygon_area);
    }

    // Rescale map to fit into a rectangular box [0, lx] * [0, ly]
    rescale_map(
        args_.n_grid_rows_or_cols,
        args_.world);

    // Store original coordinates for morphing animation
    if (args_.redirect_exports_to_stdout) {
        store_original_geo_divs();
    }

    if (args_.simplify) {
        std::cerr << "Start of initial simplification of " << pos_
                << std::endl;

        // Simplification reduces the number of points used to represent the
        // GeoDivs in the inset, thereby reducing output file sizes and
        // run-times
        simplify(args_.target_points_per_inset);
        std::cerr << "End of initial simplification of " << pos_ << std::endl;
    }

    // Plot if requested
    if (args_.plot_polygons) {

      // Color if necessary
      auto_color();
      write_cairo_map(
          inset_name_ + "_input",
          args_.plot_grid);
    }
}