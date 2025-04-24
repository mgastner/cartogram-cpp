#include "inset_state.hpp"
#include "write_image.hpp"

// ======================== Image Coloring ========================

// Functions deal with colors
Color interpolate_color(
  const double x,
  const double xmin,
  const double xmax,
  const Color ymin,
  const Color ymax)
{
  Color interpolated_color;
  // Interpolate color for red, green and blue value
  for (char c : {'r', 'g', 'b'}) {
    interpolated_color(c) =
      ((x - xmin) * ymax(c) + (xmax - x) * ymin(c)) / (xmax - xmin);
  }
  return interpolated_color;
}

// Diverging colour palette, mean accounted for
Color heatmap_color(
  const double dens,
  const double dens_min,
  const double dens_mean,
  const double dens_max)
{

  // Assign possible categories for red, green, blue
  const std::vector<Color> colors = {

    // Red to blue
    // Color("#67001f"),
    // Color("#b2182b"),
    // Color("#d6604d"),
    // Color("#f4a582"),
    // Color("#fddbc7"),
    // Color("#ffffff"),
    // Color("#d1e5f0"),
    // Color("#92c5de"),
    // Color("#4393c3"),
    // Color("#2166ac"),
    // Color("#053061")

    // Red to yellow to blue
    // Color("#a50026"),
    // Color("#d73027"),
    // Color("#f46d43"),
    // Color("#fdae61"),
    // Color("#fee090"),
    // Color("#ffffbf"),
    // Color("#e0f3f8"),
    // Color("#abd9e9"),
    // Color("#74add1"),
    // Color("#4575b4"),
    // Color("#313695")

    // Turqoise to brown
    Color("#543005"),
    Color("#8c510a"),
    Color("#bf812d"),
    Color("#dfc27d"),
    Color("#f6e8c3"),

    // Mean density
    Color("#ffffff"),  // white
    // Color("#f5f5f5"), // off-white

    Color("#c7eae5"),
    Color("#80cdc1"),
    Color("#35978f"),
    Color("#01665e"),
    Color("#003c30")

    // Original
    // Color(0.33, 0.19, 0.02),
    // Color(0.55, 0.32, 0.04),
    // Color(0.75, 0.51, 0.18),
    // Color(0.87, 0.76, 0.49),
    // Color(0.96, 0.91, 0.76),
    // Color(0.99, 0.96, 0.89),
    // Color(0.78, 0.92, 0.90),
    // Color(0.50, 0.80, 0.76),
    // Color(0.21, 0.59, 0.56),
    // Color(0.00, 0.40, 0.37),
    // Color(0.00, 0.24, 0.19)
  };
  int n_categories = colors.size();
  double xmin, xmax;
  int color_category;

  // If no discernible difference between dens and miniimum density, set
  // lowest
  if (std::fabs(dens - dens_min) <= dbl_resolution) {
    return colors[n_categories - 1];
  }

  // Choose color category
  if (dens >= dens_max) {
    return colors[0];
  } else if (dens > dens_mean) {
    color_category = 5 * (dens_max - dens) / (dens_max - dens_mean);
    xmax = dens_max - 0.2 * color_category * (dens_max - dens_mean);
    xmin = xmax - 0.2 * (dens_max - dens_mean);

    // Assign color category 0 if dens_max and dens are very close
    color_category = std::max(color_category, 0);
  } else if (dens > dens_min) {
    color_category = 5 * (dens_mean - dens) / (dens_mean - dens_min) + 5;
    xmax = dens_mean - 0.2 * (color_category - 5) * (dens_mean - dens_min);
    xmin = xmax - 0.2 * (dens_mean - dens_min);

    // Assign color category 9 if dens_min and dens are very close
    color_category = std::min(color_category, n_categories - 2);
  } else {
    return colors[n_categories - 1];
  }
  return interpolate_color(
    dens,
    xmin,
    xmax,
    colors[color_category + 1],
    colors[color_category]);
}

// Sequential colour palette, mean not accounted for
Color grid_cell_color(
  const double area,
  const double max_area,
  const double min_area)
{
  // Assign possible categories for red, green, blue
  const std::vector<Color> colors = {

    // White to purple
    Color("#ffffff"),
    Color("#efedf5"),
    Color("#dadaeb"),
    Color("#bcbddc"),
    Color("#9e9ac8"),
    Color("#807dba"),
    Color("#6a51a3"),
    Color("#54278f"),
    Color("#3f007d")

    // White to green
    // Color("#ffffff"),
    // Color("#e5f5e0"),
    // Color("#c7e9c0"),
    // Color("#a1d99b"),
    // Color("#74c476"),
    // Color("#41ab5d"),
    // Color("#238b45"),
    // Color("#006d2c"),
    // Color("#00441b")

    // // White to red
    // Color(1.000, 0.961, 0.941),
    // Color(0.996, 0.878, 0.824),
    // Color(0.988, 0.733, 0.631),
    // Color(0.988, 0.572, 0.447),
    // Color(0.984, 0.416, 0.290),
    // Color(0.937, 0.231, 0.173),
    // Color(0.796, 0.094, 0.114),
    // Color(0.647, 0.058, 0.082),
    // Color(0.404, 0.000, 0.050)
  };
  int n_categories = colors.size();

  // Normalize area to [0,1] and make it logarithmic
  double ratio = (log(area) - log(min_area)) / (log(max_area) - log(min_area));

  // Determine color category
  double category = fmax(0, ratio * (n_categories - 1) - 10e-6);
  int xmin = floor(category);
  int xmax = std::ceil(category);

  if (area >= max_area) {
    return colors[n_categories - 1];
  } else if (area == min_area) {
    return colors[0];
  } else {
    Color interpolated_color;
    for (char c : {'r', 'g', 'b'}) {
      interpolated_color(c) =
        colors[xmin](c) +
        (colors[xmax](c) - colors[xmin](c)) * (category - xmin);
    }
    return interpolated_color;
  }
}

// Calculates all grid cell colors from cartogram map to be used in
// the equal area grid heatmap
std::vector<std::vector<Color>> grid_cell_colors(
  unsigned int cell_width,
  InsetState &inset_state)
{
  // Initialize max and min area
  double max_area, min_area;
  std::tie(max_area, min_area) =
    max_and_min_grid_cell_area(cell_width, inset_state);

  // Initialize colors
  std::vector<std::vector<Color>> colors(
    inset_state.lx() - cell_width,
    std::vector<Color>(inset_state.ly() - cell_width));

  // Iterate over grid cells
  for (unsigned int i = 0; i < inset_state.lx() - cell_width;
       i += cell_width) {
    for (unsigned int j = 0; j < inset_state.ly() - cell_width;
         j += cell_width) {
      const double area = grid_cell_area(i, j, cell_width, inset_state);
      colors[i][j] = grid_cell_color(area, max_area, min_area);
    }
  }
  return colors;
}

void InsetState::write_grid_colors_on_surface(
  cairo_t *cr,
  bool plot_equal_area_map,
  bool crop_polygons)
{
  unsigned int cell_width = 1;

  // Get colors
  const auto colors = grid_cell_colors(cell_width, *this);

  // Set line width of grid lines
  cairo_set_line_width(cr, 5e-6 * std::min(lx_, ly_));

  // Print the color of all grid cells and exit
  if (!crop_polygons) {

    // Iterate over grid cells
    for (unsigned int i = 0; i < lx_ - cell_width; i += cell_width) {
      for (unsigned int j = 0; j < ly_ - cell_width; j += cell_width) {
        // Set color of the border of the grid polygon
        cairo_set_source_rgb(
          cr,
          colors[i][j].r,
          colors[i][j].g,
          colors[i][j].b);

        // Draw grid cell by connecting edge points
        const Polygon cell_edge_points =
          grid_cell_edge_points(i, j, cell_width, plot_equal_area_map);
        cairo_move_to(
          cr,
          cell_edge_points[0].x(),
          ly_ - cell_edge_points[0].y());
        for (unsigned int k = 1; k < cell_edge_points.size(); ++k) {
          cairo_line_to(
            cr,
            cell_edge_points[k].x(),
            ly_ - cell_edge_points[k].y());
        }

        // Fill the grid polygon with color
        cairo_fill_preserve(cr);
        cairo_stroke(cr);
      }
    }

    // Don't plot uncropped version
    return;
  }

  // Draw cartogram polygons or equal area map polygons
  const std::vector<GeoDiv> &geo_divs =
    plot_equal_area_map ? geo_divs_original_ : geo_divs_;

  const std::vector<GeoDiv> &geo_divs_original = geo_divs_original_;
  std::map<std::string, Bbox> original_bboxes;

  for (const auto &gd : geo_divs_original) {
    original_bboxes[gd.id()] = gd.bbox();
  }

  // Clip to shape
  for (const auto &gd : geo_divs) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto ext_ring = pwh.outer_boundary();
      cairo_new_path(cr);
      cairo_move_to(cr, ext_ring[0].x(), ly_ - ext_ring[0].y());

      // Plot each point in exterior ring
      for (unsigned int i = 1; i < ext_ring.size(); ++i) {
        cairo_line_to(cr, ext_ring[i].x(), ly_ - ext_ring[i].y());
      }

      // Close entire path
      cairo_close_path(cr);
      cairo_clip(cr);

      Bbox gd_bbox = original_bboxes.at(gd.id());

      // Iterate over grid cells
      for (unsigned int i = gd_bbox.xmin() - 1; i < gd_bbox.xmax() + 1;
           i += cell_width) {
        for (unsigned int j = gd_bbox.ymin() - 1; j < gd_bbox.ymax() + 1;
             j += cell_width) {

          // Set color of the border of the grid polygon
          cairo_set_source_rgb(
            cr,
            colors[i][j].r,
            colors[i][j].g,
            colors[i][j].b);

          // Draw grid cells by connecting edge points
          const auto cell_edge_points =
            grid_cell_edge_points(i, j, cell_width, plot_equal_area_map);

          // Move to first grid edge
          cairo_move_to(
            cr,
            cell_edge_points[0].x(),
            ly_ - cell_edge_points[0].y());

          // Draw remaining grid cells
          for (unsigned int k = 1; k < cell_edge_points.size(); ++k) {
            cairo_line_to(
              cr,
              cell_edge_points[k].x(),
              ly_ - cell_edge_points[k].y());
          }

          // Fill the grid polygon with color
          cairo_fill_preserve(cr);
          cairo_stroke(cr);
        }
      }

      // Remove GeoDiv clip
      cairo_reset_clip(cr);
    }
  }
}