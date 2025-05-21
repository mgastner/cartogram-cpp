#ifndef USE_CAIRO
#define USE_CAIRO 0
#endif

#if USE_CAIRO

void right_aligned_text(
  cairo_t *cr,
  double number,
  double x,
  double y,
  double font_size)
{
  std::string text = std::to_string(static_cast<int>(round(number)));
  cairo_text_extents_t extents;
  cairo_text_extents(cr, text.c_str(), &extents);
  double move = extents.width;
  cairo_move_to(cr, x - move - font_size * 0.75, y);
  cairo_show_text(cr, text.c_str());
}

// Returns grid cell area based on edge points
static double grid_cell_area(
  unsigned int x,
  unsigned int y,
  unsigned int cell_width,
  InsetState &inset_state)
{
  // Taking absolule to ensure we get the area irrespective of direction
  return abs(inset_state.grid_cell_edge_points(x, y, cell_width).area());
}

void InsetState::write_grid_heatmap_data(const std::string filename)
{
  unsigned int cell_width = 1;
  const unsigned int resolution = 1;

  InsetState is_copy = (*this);
  is_copy.set_geo_divs(geo_divs_original_);
  auto intersections_with_rays =
    is_copy.intersec_with_parallel_to('x', resolution);
  std::vector<std::vector<double>> exists(lx_, std::vector<double>(ly_, 0));

  // Mark all squares that are inside map with 1
  for (unsigned int y = 0; y < ly_; y += 1.0) {

    // Intersections for one ray
    auto intersections_at_y = intersections_with_rays[y];

    // Sort intersections in ascending order
    std::sort(intersections_at_y.begin(), intersections_at_y.end());

    // Fill GeoDivs by iterating over intersections
    for (unsigned int i = 0; i < intersections_at_y.size(); i += 2) {
      const double left_x = intersections_at_y[i].x();
      const double right_x = intersections_at_y[i + 1].x();

      // Fill each cell between intersections
      for (unsigned int m = ceil(left_x); m <= ceil(right_x); ++m) {
        exists[m - 1][y] += 1;
      }
    }
  }

  std::ofstream f_csv;
  f_csv.open(filename);

  // Fill rho_init with the ratio of rho_num to exists
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      if (exists[i][j]) {
        const double area = grid_cell_area(i, j, cell_width, *this);
        f_csv << i << ", " << j << ", " << area << "\n";
      }
    }
  }
  std::cerr << "Grid heatmap data written " + filename << std::endl;
}

void write_grid_heatmap_bar_on_surface(
  double min_value,
  double max_value,
  cairo_t *cr,
  Bbox bbox_bar,
  std::vector<std::pair<double, double>> major_ticks,
  std::vector<std::pair<double, double>> minor_ticks,
  const unsigned int ly)
{
  const int n_gradient_bars = 500;

  // get bar coordinates
  const double xmin_bar = bbox_bar.xmin();
  const double xmax_bar = bbox_bar.xmax();
  const double ymin_bar = bbox_bar.ymin();
  const double ymax_bar = bbox_bar.ymax();

  const double bar_width = xmax_bar - xmin_bar;
  const double bar_height = ymax_bar - ymin_bar;

  // Draw the outer bar lines
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_set_line_width(cr, 1.0);
  cairo_rectangle(cr, xmin_bar, ymin_bar, bar_width, -bar_height);
  cairo_stroke(cr);

  // Calculate individual bar gradient segment property
  const double gradient_segment_height =
    (ymax_bar - ymin_bar) / n_gradient_bars;
  const double gradient_segment_value =
    (max_value - min_value) / n_gradient_bars;

  // Draw the gradient segment rectangles
  double value_at_gradient_segment = min_value;
  double overlap = 0.1;

  for (double y = ymin_bar; y <= ymax_bar; y += gradient_segment_height) {
    Color color = grid_cell_color(
      exp(value_at_gradient_segment),
      exp(max_value),
      exp(min_value));
    cairo_set_source_rgb(cr, color.r, color.g, color.b);
    cairo_rectangle(
      cr,
      xmin_bar,
      ly - y - overlap,
      bar_width,
      gradient_segment_height + overlap);
    cairo_fill(cr);
    value_at_gradient_segment += gradient_segment_value;
  }

  unsigned int font_size = 10;

  // Set font properties
  cairo_select_font_face(
    cr,
    "Sans",
    CAIRO_FONT_SLANT_NORMAL,
    CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size(cr, font_size);

  // Tick width outside bar
  double half_tick_width = 4.0;

  // Draw the ticks and nice_numbers
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_set_line_width(cr, 1.1);
  for (auto tick : major_ticks) {
    double area = tick.first;
    int NiceNumber = tick.second;
    if (area > min_value and area < max_value) {
      double y =
        ((log(area) - log(min_value)) / (log(max_value) - log(min_value))) *
          (ymax_bar - ymin_bar) +
        ymin_bar;
      cairo_move_to(cr, xmin_bar + half_tick_width, ly - y);
      cairo_line_to(cr, xmin_bar - half_tick_width, ly - y);

      // Right-align text
      right_aligned_text(
        cr,
        NiceNumber,
        xmin_bar - half_tick_width,
        ly - y + font_size / 2.0,
        font_size);
    }
  }
  cairo_set_line_width(cr, 0.75);
  half_tick_width /= 2.0;
  for (auto ticks : minor_ticks) {
    double area = ticks.first;
    if (area > min_value and area < max_value) {
      double y =
        ((log(area) - log(min_value)) / (log(max_value) - log(min_value))) *
          (ymax_bar - ymin_bar) +
        ymin_bar;

      cairo_move_to(cr, xmin_bar + half_tick_width, ly - y);
      cairo_line_to(cr, xmin_bar - half_tick_width, ly - y);
      cairo_stroke(cr);
    }
  }
}

double grid_cell_target_area_per_km(
  const unsigned int i,
  const unsigned int j,
  const double total_target_area,
  const double total_inset_area,
  InsetState &inset_state)
{
  const double cell_target_area = grid_cell_target_area(
    i,
    j,
    total_target_area,
    total_inset_area,
    inset_state);
  const double cell_area_km = grid_cell_area_km(inset_state, i, j);

  const double cell_target_area_per_km = cell_target_area / cell_area_km;

  return cell_target_area_per_km;
}

// Returns the index of largest and smallest grid cell area to be used for
// grid heatmap generation
static std::pair<Point, Point> max_and_min_grid_cell_area_index(
  unsigned int cell_width,
  InsetState &inset_state)
{
  // Initialize max and min area
  double max_area = -dbl_inf;
  double min_area = dbl_inf;
  unsigned int max_i = 0;
  unsigned int max_j = 0;
  unsigned int min_i = 0;
  unsigned int min_j = 0;

  // Iterate over grid cells
  for (unsigned int i = 0; i < inset_state.lx() - cell_width;
       i += cell_width) {
    for (unsigned int j = 0; j < inset_state.ly() - cell_width;
         j += cell_width) {
      const double area = grid_cell_area(i, j, cell_width, inset_state);
      if (area > max_area) {
        max_area = area;
        max_i = i;
        max_j = j;
      }
      if (area < min_area) {
        min_area = area;
        min_i = i;
        min_j = j;
      }
    }
  }
  return std::make_pair(Point(max_i, max_j), Point(min_i, min_j));
}

// Outputs a SVG/PS file of grid heatmap
void write_grid_heatmap_image(
  const std::string filename,
  const bool plot_equal_area_map,
  const bool crop_polygons,
  InsetState &inset_state)
{
  // Whether to draw bar on the cairo surface
  const bool draw_bar = false;

  // Create a cairo surface
  cairo_surface_t *surface;
  surface = cairo_svg_surface_create(
    filename.c_str(),
    inset_state.lx(),
    inset_state.ly());
  cairo_t *cr = cairo_create(surface);
  add_white_background(cr);

  // Get inset areas
  const double total_ta = inset_state.total_target_area();
  const double total_ia = inset_state.total_inset_area();

  // const Bbox bbox_bar = get_bbox_bar(15, 150);

  // Get the max and min grid cell area points
  Point max_area_cell_point, min_area_cell_point;

  std::tie(max_area_cell_point, min_area_cell_point) =
    max_and_min_grid_cell_area_index(1, inset_state);

  const double max_area_cell_point_area = grid_cell_area(
    max_area_cell_point.x(),
    max_area_cell_point.y(),
    1,
    inset_state);
  const double min_area_cell_point_area = grid_cell_area(
    min_area_cell_point.x(),
    min_area_cell_point.y(),
    1,
    inset_state);

  // Get the max and min grid cell target area per km
  double max_target_area_per_km = grid_cell_target_area_per_km(
    max_area_cell_point.x(),
    max_area_cell_point.y(),
    total_ta,
    total_ia,
    inset_state);
  double min_target_area_per_km = grid_cell_target_area_per_km(
    min_area_cell_point.x(),
    min_area_cell_point.y(),
    total_ta,
    total_ia,
    inset_state);

  if (plot_equal_area_map) {
    std::cerr << std::endl;
    std::cerr << "Max target area per km: " << max_target_area_per_km
              << std::endl;
    std::cerr << "Min target area per km: " << min_target_area_per_km
              << std::endl
              << std::endl;
  }

  std::vector<int> nice_numbers =
    get_nice_numbers_for_bar(max_target_area_per_km);

  std::vector<std::pair<double, double>> major_ticks, minor_ticks;

  std::tie(major_ticks, minor_ticks) = get_ticks(
    9,
    min_target_area_per_km,
    max_target_area_per_km,
    min_area_cell_point_area,
    max_area_cell_point_area,
    nice_numbers);

  // Draw colors
  inset_state.write_grid_colors_on_surface(
    cr,
    plot_equal_area_map,
    crop_polygons);

  // Draw polygons without color
  write_polygons_on_surface(cr, false, false, inset_state);

  if (plot_equal_area_map)
    write_grid_on_surface(cr, inset_state);

  // trim_grid_heatmap(cr, 20);

  if (draw_bar) {
    std::string bar_filename = "bar_" + filename;

    // Create a cairo bar_surface
    cairo_surface_t *bar_surface =
      cairo_svg_surface_create(bar_filename.c_str(), 160, 400);
    cairo_t *bar_cr = cairo_create(bar_surface);

    // Write bar
    write_grid_heatmap_bar_on_surface(
      min_area_cell_point_area,
      max_area_cell_point_area,
      bar_cr,
      Bbox(75, 200, 95, 350),
      major_ticks,
      minor_ticks,
      400);

    cairo_show_page(bar_cr);
    cairo_surface_destroy(bar_surface);
    cairo_destroy(bar_cr);
  }
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

// Function to show the density bar on the cairo surface
void write_density_bar_on_surface(
  const double min_value,
  const double mean_value,
  const double max_value,
  cairo_t *cr,
  Bbox bbox_bar,
  const unsigned int ly)
{
  const int n_gradient_bars = 500;

  // get bar coordinates
  const double xmin_bar = bbox_bar.xmin();
  const double xmax_bar = bbox_bar.xmax();
  const double ymin_bar = bbox_bar.ymin();
  const double ymax_bar = bbox_bar.ymax();

  const double bar_width = xmax_bar - xmin_bar;
  const double bar_height = ymax_bar - ymin_bar;

  // Draw the outer bar lines
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_set_line_width(cr, 1.0);
  cairo_rectangle(cr, xmin_bar, ymin_bar, bar_width, -bar_height);
  cairo_stroke(cr);

  // position of mean line along bar
  // const double ymean_bar = ((ymax_bar - ymin_bar) / (max_value - min_value))
  // *
  //                            (mean_value - min_value) +
  //                          ymin_bar;

  // calculate individual bar gradient segment property
  const double gradient_segment_height =
    (ymax_bar - ymin_bar) / n_gradient_bars;
  const double gradient_segment_value =
    abs(max_value - min_value) / n_gradient_bars;

  // Draw the gradient segment rectangles
  double value_at_gradient_segment = min_value;
  double overlap = 0.1;

  for (double y = ymin_bar; y <= ymax_bar; y += gradient_segment_height) {
    Color color = heatmap_color(
      value_at_gradient_segment,
      min_value,
      mean_value,
      max_value);
    cairo_set_source_rgb(cr, color.r, color.g, color.b);
    cairo_rectangle(
      cr,
      xmin_bar,
      ly - y - overlap,
      bar_width,
      gradient_segment_height + overlap);
    cairo_fill(cr);
    value_at_gradient_segment += gradient_segment_value;
  }

  // Draw the mean line
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  // cairo_set_line_width(cr, 0.9);
  // cairo_move_to(cr, xmax_bar - bar_width / 2, ly - ymean_bar);
  // cairo_line_to(cr, xmin_bar, ly - ymean_bar);
  // cairo_stroke(cr);

  unsigned int font_size = 10;

  // Set font properties
  cairo_select_font_face(
    cr,
    "Sans",
    CAIRO_FONT_SLANT_NORMAL,
    CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size(cr, font_size);

  // Write residual density
  // cairo_move_to(
  //   cr,
  //   xmin_bar - bar_width / 2 - 1,
  //   ly - ymax_bar - (font_size * 3.0));
  // cairo_show_text(cr, "Residual");
  // cairo_move_to(
  //   cr,
  //   xmin_bar - bar_width / 2 - 1,
  //   ly - ymax_bar - (font_size * 2.0));
  // cairo_show_text(cr, "Density (km-²)");

  // Tick width outside bar
  double half_tick_width = bar_width / 8.0;

  // Write max value at top of bar
  right_aligned_text(
    cr,
    std::floor(max_value),
    xmin_bar - half_tick_width,
    ly - ymax_bar,
    font_size);

  // Write mean value
  // right_aligned_text(
  //   cr,
  //   mean_value,
  //   xmin_bar,
  //   ly - ymean_bar + (font_size / 4.0),
  //   font_size);

  // Write min_value
  right_aligned_text(
    cr,
    std::ceil(min_value),
    xmin_bar - half_tick_width,
    ly - ymin_bar + (font_size * 1.2),
    font_size);

  long long magnitude =
    std::max(std::pow(10.0, std::floor(std::log10(max_value))), 1.0);

  long long original_mag = magnitude;

  if (max_value > magnitude * 4 || min_value < -magnitude * 4) {
    magnitude *= 2;
  } else if (
    max_value < magnitude * 2 && std::abs(min_value) < magnitude * 2) {
    magnitude /= 2;
  }

  // Highest magnitude
  long long highest_mag = magnitude * (std::floor(max_value / magnitude));

  if (highest_mag > max_value - (0.4 * original_mag)) {
    highest_mag -= magnitude;
  }

  if (magnitude >= 0.1) {
    cairo_set_line_width(cr, 0.7);
    double bar_ratio = (ymax_bar - ymin_bar) / (max_value - min_value);

    // Ticks
    for (long long i = highest_mag;
         i - (0.5 * magnitude) > min_value || i == mean_value;
         i -= magnitude) {
      double tick = bar_ratio * (i - min_value) + ymin_bar;
      // if (i == mean_value) {
      //   tick += (font_size / 4.0);
      // }
      cairo_move_to(cr, xmin_bar + half_tick_width, ly - tick);
      cairo_line_to(cr, xmin_bar - half_tick_width, ly - tick);
      cairo_stroke(cr);
      right_aligned_text(
        cr,
        i,
        xmin_bar - half_tick_width,
        ly - tick + (font_size / 2.0),
        font_size);
    }
  }
}

void InsetState::write_density_image(
  const std::string filename,
  const bool plot_pycnophylactic)
{
  double *density = rho_init_.as_1d_array();
  std::cerr << "Writing " << filename << std::endl;
  // Whether to draw bar on the cairo surface
  const bool draw_bar = false;
  cairo_surface_t *surface =
    cairo_svg_surface_create(filename.c_str(), lx_, ly_);

  // Create a cairo surface
  cairo_t *cr = cairo_create(surface);
  cairo_set_line_width(cr, 0);

  // Determine range of densities
  const double dens_min = dens_min_;
  const double dens_mean = dens_mean_;
  const double dens_max = dens_max_;
  const double exterior_density = exterior_density_;
  const double each_grid_cell_area_km = grid_cell_area_km(*this, 0, 0);

  // Debugging prints
  std::cerr << "dens_min: " << dens_min << std::endl;
  std::cerr << "dens_mean: " << dens_mean << std::endl;
  std::cerr << "dens_max: " << dens_max << std::endl;
  std::cerr << "exterior_density: " << exterior_density << std::endl;
  std::cerr << "each_grid_cell_area_km: " << each_grid_cell_area_km
            << std::endl;

  // Crop it too
  if (plot_pycnophylactic) {
    unsigned int cell_width = 1;
    add_white_background(cr);

    // Clip to shape and print in sequential scale
    for (const auto &gd : geo_divs_) {
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

        Bbox gd_bbox = gd.bbox();

        // Iterate over grid cells
        for (unsigned int i = gd_bbox.xmin(); i < gd_bbox.xmax();
             i += cell_width) {
          for (unsigned int j = gd_bbox.ymin(); j < gd_bbox.ymax();
               j += cell_width) {

            double target_area_km =
              density[i * ly_ + j] / each_grid_cell_area_km;

            // Values here used are "Max target area per km" and
            // "Min target area per km", which is obtained by running the
            // code with the "plot_pycnophylactic" -y flag set to true
            // Update here for new map
            Color color = grid_cell_color(target_area_km, 660.058, 0.660816);

            // Get four points of the square
            double x_min = i - 0.5 * sq_overlap;
            double y_min = j - 0.5 * sq_overlap;
            double x_max = i + 1 + 0.5 * sq_overlap;
            double y_max = j + 1 + 0.5 * sq_overlap;

            cairo_move_to(cr, x_min, ly_ - y_min);
            cairo_line_to(cr, x_max, ly_ - y_min);
            cairo_line_to(cr, x_max, ly_ - y_max);
            cairo_line_to(cr, x_min, ly_ - y_max);

            cairo_set_source_rgb(cr, color.r, color.g, color.b);
            cairo_fill(cr);
            cairo_set_source_rgb(cr, 0, 0, 0);
            cairo_stroke(cr);

            // Fill the grid polygon with color
            cairo_fill_preserve(cr);
            cairo_stroke(cr);
          }
        }

        // Remove GeoDiv clip
        cairo_reset_clip(cr);
      }
    }
  } else {

    // Set background color to that of exterior_density
    Color exterior_density_color =
      heatmap_color(exterior_density, dens_min, dens_mean, dens_max);
    cairo_set_source_rgb(
      cr,
      exterior_density_color.r / 255.0,
      exterior_density_color.g / 255.0,
      exterior_density_color.b / 255.0);
    cairo_paint(cr);
    for (unsigned int i = 0; i < lx_; ++i) {
      for (unsigned int j = 0; j < ly_; ++j) {

        Color color =
          heatmap_color(density[i * ly_ + j], dens_min, dens_mean, dens_max);

        // Skip plotting exterior_density color
        if (color == exterior_density_color)
          continue;

        // Get four points of the square
        double x_min = i - 0.5 * sq_overlap;
        double y_min = j - 0.5 * sq_overlap;
        double x_max = i + 1 + 0.5 * sq_overlap;
        double y_max = j + 1 + 0.5 * sq_overlap;

        cairo_move_to(cr, x_min, ly_ - y_min);
        cairo_line_to(cr, x_max, ly_ - y_min);
        cairo_line_to(cr, x_max, ly_ - y_max);
        cairo_line_to(cr, x_min, ly_ - y_max);

        cairo_set_source_rgb(
          cr,
          color.r / 255.0,
          color.g / 255.0,
          color.b / 255.0);
        cairo_fill(cr);
        cairo_set_source_rgb(cr, 0, 0, 0);
        cairo_stroke(cr);

        // Write the density value at the center of the square
        // cairo_set_source_rgb(cr, 0, 0, 0);
        // cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL,
        // CAIRO_FONT_WEIGHT_NORMAL); cairo_set_font_size(cr, 10); std::string
        // density_value = std::to_string(density[i * ly_ + j]);
        // cairo_text_extents_t extents;
        // cairo_text_extents(cr, density_value.c_str(), &extents);
        // double x_center = (x_min + x_max) / 2 - extents.width / 2;
        // double y_center = ly_ - (y_min + y_max) / 2 + extents.height / 2;
        // cairo_move_to(cr, x_center, y_center);
        // cairo_show_text(cr, density_value.c_str());
      }
    }
  }
  write_quadtree_rectangles_on_surface(
    cr,
    quadtree_bboxes_,
    Color{0.6, 0.6, 0.6},
    ly_);
  write_polygons_on_surface(cr, false, false, *this, 0.025);

  if (draw_bar && !plot_pycnophylactic) {
    std::string bar_filename = "bar_" + filename;

    // Create a cairo bar_surface
    cairo_surface_t *bar_surface =
      cairo_svg_surface_create(bar_filename.c_str(), 160, 400);

    cairo_t *bar_cr = cairo_create(bar_surface);

    write_density_bar_on_surface(
      dens_min - dens_mean,
      0,
      dens_max - dens_mean,
      bar_cr,
      Bbox(75, 200, 95, 350),
      400);

    cairo_show_page(bar_cr);
    cairo_surface_destroy(bar_surface);
    cairo_destroy(bar_cr);
  }

  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

// Given padding, makes padding area white
void trim_grid_heatmap(cairo_t *cr, double padding, InsetState &inset_state)
{
  // Canvas dimension
  Bbox is_bbox = inset_state.bbox();
  double xmin = is_bbox.xmin() - padding;
  double xmax = is_bbox.xmax() + padding;
  double ymin = is_bbox.ymin() - padding;
  double ymax = is_bbox.ymax() + padding;

  // Color white outside is_bbox
  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_rectangle(cr, 0, 0, xmin, ymax);
  cairo_rectangle(cr, xmax, 0, inset_state.lx() - xmax, ymax);
  cairo_rectangle(cr, 0, ymax, inset_state.lx(), inset_state.ly() - ymax);
  cairo_rectangle(cr, 0, 0, inset_state.lx(), ymin);
  cairo_fill(cr);
}

std::vector<int> get_nice_numbers_for_bar(const double max_target_area_per_km)
{
  std::vector<int> nice_numbers;
  int NiceNumber = 1;
  nice_numbers.push_back(NiceNumber);
  while (NiceNumber < max_target_area_per_km) {
    NiceNumber = NiceNumber * 10;
    nice_numbers.push_back(NiceNumber);
  }
  return nice_numbers;
}

std::vector<std::pair<double, double>> get_major_ticks(
  const double min_target_area_per_km,
  const double max_target_area_per_km,
  const double min_area_cell_point_area,
  const double max_area_cell_point_area,
  std::vector<int> nice_numbers)

{
  std::vector<std::pair<double, double>> ticks;
  for (auto niceNumber : nice_numbers) {
    double NiceNumberRatio = (niceNumber - min_target_area_per_km) /
                             (max_target_area_per_km - min_target_area_per_km);
    double area =
      min_area_cell_point_area +
      NiceNumberRatio * (max_area_cell_point_area - min_area_cell_point_area);
    ticks.push_back(std::make_pair(area, niceNumber));
  }
  return ticks;
}

std::vector<std::pair<double, double>> get_minor_ticks(
  int n_ticks_per_major,
  const double min_target_area_per_km,
  const double max_target_area_per_km,
  const double min_area_cell_point_area,
  const double max_area_cell_point_area,
  std::vector<int> nice_numbers)

{
  n_ticks_per_major += 2;
  std::vector<std::pair<double, double>> minor_ticks;
  const int n_major_ticks = nice_numbers.size();
  for (int i = 0; i < n_major_ticks - 1; i++) {
    double first_major_tick = (double)nice_numbers[i];
    double second_major_tick = (double)nice_numbers[i + 1];
    double minor_tick_ratio;
    if (first_major_tick == 1.0) {
      minor_tick_ratio =
        (second_major_tick - first_major_tick) / (n_ticks_per_major - 2);
    } else {
      minor_tick_ratio =
        (second_major_tick - first_major_tick) / (n_ticks_per_major - 1);
    }
    for (int j = 1; j < n_ticks_per_major - 1; j++) {
      double minor_tick = first_major_tick + minor_tick_ratio * j;
      double NiceNumberRatio =
        (minor_tick - min_target_area_per_km) /
        (max_target_area_per_km - min_target_area_per_km);
      double area = min_area_cell_point_area +
                    NiceNumberRatio *
                      (max_area_cell_point_area - min_area_cell_point_area);
      minor_ticks.push_back(std::make_pair(area, minor_tick));
    }
  }
  return minor_ticks;
}

std::pair<
  std::vector<std::pair<double, double>>,
  std::vector<std::pair<double, double>>>
get_ticks(
  const int n_ticks_per_major,
  const double min_target_area_per_km,
  const double max_target_area_per_km,
  const double min_area_cell_point_area,
  const double max_area_cell_point_area,
  std::vector<int> nice_numbers)
{
  std::vector<std::pair<double, double>> major_ticks = get_major_ticks(
    min_target_area_per_km,
    max_target_area_per_km,
    min_area_cell_point_area,
    max_area_cell_point_area,
    nice_numbers);
  std::vector<std::pair<double, double>> minor_ticks = get_minor_ticks(
    n_ticks_per_major,
    min_target_area_per_km,
    max_target_area_per_km,
    min_area_cell_point_area,
    max_area_cell_point_area,
    nice_numbers);
  std::pair<
    std::vector<std::pair<double, double>>,
    std::vector<std::pair<double, double>>>
    ticks(major_ticks, minor_ticks);
  return ticks;
}

Bbox get_bbox_bar(
  const double bar_width,
  const double bar_height,
  InsetState &inset_state)
{
  const Bbox is_bbox = inset_state.bbox();

  // Position the bar 25 pixels to the right of the is_bbox
  const double xmin_bar = is_bbox.xmax() + 35;
  const double xmax_bar = xmin_bar + bar_width;

  // Position the bar at the middle of the is_bbox y coordinates
  double ymid_bar = (is_bbox.ymax() + is_bbox.ymin()) / 2 - 25;
  double ymin_bar = ymid_bar - bar_height / 2;
  double ymax_bar = ymid_bar + bar_height / 2;

  return Bbox(xmin_bar, ymin_bar, xmax_bar, ymax_bar);
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

      Bbox gd_bbox;
      try {
        gd_bbox = original_bboxes.at(gd.id());
      } catch (const std::out_of_range &e) {
        std::cerr << "ERROR: Key '" << gd.id()
                  << "' not found in original_bboxes. "
                  << "Exception: " << e.what() << std::endl;
        // Re-throw, or return a default value
        throw;
      }

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

#endif
