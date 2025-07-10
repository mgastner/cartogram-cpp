#pragma once

#include "colors.hpp"
#include "inset_state.hpp"

// ======================== Image Coloring ========================

// Functions deal with colors
static Color interpolate_color(
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
static Color heatmap_color(
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
  const size_t n_categories = colors.size();
  assert(n_categories > 0 && "Color categories must be defined");
  double xmin, xmax;
  size_t color_category;

  // If no discernible difference between dens and miniimum density, set
  // lowest
  if (std::fabs(dens - dens_min) <= dbl_resolution) {
    return colors[n_categories - 1];
  }

  // Choose color category
  if (dens >= dens_max) {
    return colors[0];
  } else if (dens > dens_mean) {
    color_category = static_cast<size_t>(
      std::max(0.0, 5 * (dens_max - dens) / (dens_max - dens_mean)));
    xmax = dens_max -
           0.2 * static_cast<double>(color_category) * (dens_max - dens_mean);
    xmin = xmax - 0.2 * (dens_max - dens_mean);
  } else if (dens > dens_min) {
    color_category = static_cast<size_t>(
      std::max(0.0, 5 * (dens_mean - dens) / (dens_mean - dens_min) + 5));
    xmax = dens_mean - 0.2 * static_cast<double>((color_category - 5)) *
                         (dens_mean - dens_min);
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
