#include "write_image.h"
#include <iostream>

void write_triangles_on_cairo_surface(cairo_t *cr, Delaunay &dt, color clr)
{
  // Draw the triangles
  for (Delaunay::Finite_faces_iterator fit = dt.finite_faces_begin();
       fit != dt.finite_faces_end();
       ++fit) {
    Point p1 = fit->vertex(0)->point();
    Point p2 = fit->vertex(1)->point();
    Point p3 = fit->vertex(2)->point();

    // set width of line
    cairo_set_line_width(cr, 0.05);

    // set color
    cairo_set_source_rgb(cr, clr.r, clr.g, clr.b);

    cairo_move_to(cr, p1.x(), 512 - p1.y());
    cairo_line_to(cr, p2.x(), 512 - p2.y());
    cairo_line_to(cr, p3.x(), 512 - p3.y());
    cairo_line_to(cr, p1.x(), 512 - p1.y());
    cairo_stroke(cr);
  }
}

void write_ps_header(const std::string filename, cairo_surface_t *surface)
{
  const std::string title = "%%Title: " + filename;
  cairo_ps_surface_dsc_comment(surface, title.c_str());
  cairo_ps_surface_dsc_comment(
    surface,
    "%%Creator: Michael T. Gastner et al.");
  cairo_ps_surface_dsc_comment(surface, "%%For: Humanity");
  cairo_ps_surface_dsc_comment(surface, "%%Copyright: License CC BY");
  cairo_ps_surface_dsc_comment(surface, "%%Magnification: 1.0000");
}

void write_point_on_cairo_surface(cairo_t *cr, Point pt, color clr)
{
  cairo_set_source_rgb(cr, clr.r, clr.g, clr.b);
  cairo_move_to(cr, pt.x(), 512 - pt.y());
  cairo_line_to(cr, pt.x() + 0.05, 512 - pt.y() - 0.05);
  cairo_close_path(cr);
  cairo_stroke(cr);
}

void InsetState::write_polygon_points_on_cairo_surface(cairo_t *cr, color clr)
{
  cairo_set_source_rgb(cr, clr.r, clr.g, clr.b);
  cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
  cairo_set_line_width(cr, 0.5);
  // Draw the shapes
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto ext_ring = pwh.outer_boundary();

      // Plot each point in exterior ring
      for (auto i : ext_ring) {
        write_point_on_cairo_surface(cr, i, clr);
      }

      // Plot holes
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
        Polygon hole = *hci;
        for (unsigned int i = 1; i <= hole.size(); ++i) {
          write_point_on_cairo_surface(cr, hole[i], clr);
        }
      }
    }
  }
}

void InsetState::write_quadtree(const std::string &filename)
{
  cairo_surface_t *surface =
    cairo_ps_surface_create((filename + ".ps").c_str(), lx_, ly_);
  cairo_t *cr = cairo_create(surface);
  write_ps_header((filename + ".ps"), surface);
  write_triangles_on_cairo_surface(cr, proj_qd_.dt, color{0.0, 0.0, 0.0});
  write_polygon_points_on_cairo_surface(cr, color{1.0, 0.0, 0.0});
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

// TODO: IS THERE A CGAL WAY OF DETERMINING WHETHER THE LABEL'S BOUNDING
//       BOX IS COMPLETELY CONTAINED IN THE POLYGON?

// TODO: SHOULD THE CRITERION FOR PRINTING A LABEL BE THAT IT FITS INSIDE THE
//       POLYGON WITH HOLES? THAT CRITERION WOULD BE MORE RESTRICTIVE THAN
//       FITTING INSIDE THE EXTERIOR RING.
bool all_points_inside_exterior_ring(
  const std::vector<Point> &pts,
  const Polygon_with_holes &pwh)
{
  for (const auto &pt : pts) {
    if (pwh.outer_boundary().has_on_unbounded_side(pt)) {
      return false;
    }
  }
  return true;
}

double font_size(
  cairo_t *cr,
  const char *label,
  const Point label_pt,
  const GeoDiv &gd)
{
  for (double fsize = max_font_size; fsize >= min_font_size; fsize -= 0.5) {
    cairo_set_font_size(cr, fsize);
    cairo_text_extents_t extents;
    cairo_text_extents(cr, label, &extents);
    const auto largest_pwh = gd.largest_polygon_with_holes();

    // Bounding box of the label
    const CGAL::Bbox_2 bb(
      label_pt.x() - 0.5 * extents.width,
      label_pt.y() - 0.5 * extents.height,
      label_pt.x() + 0.5 * extents.width,
      label_pt.y() + 0.5 * extents.height);

    // Vector of bounding-box edge points
    std::vector<Point> bb_edge_points;
    for (unsigned int i = 0; i <= 1; ++i) {
      for (unsigned int j = 0; j <= 5; ++j) {
        bb_edge_points.emplace_back(
          (j * bb.xmin() + (5 - j) * bb.xmax()) / 5,
          (i * bb.ymin() + (1 - i) * bb.ymax()));
      }
    }
    if (all_points_inside_exterior_ring(bb_edge_points, largest_pwh)) {
      return fsize;
    }
  }
  return 0.0;
}


void InsetState::write_polygons_to_cairo_surface(
  cairo_t *cr,
  const bool fill_polygons,
  const bool colors,
  const bool plot_equal_area_map)
{
  cairo_set_line_width(cr, 1e-3 * std::min(lx_, ly_));

  // Draw cartogram polygons or equal area map polygons
  const std::vector<GeoDiv> &geo_divs =
    plot_equal_area_map ? geo_divs_original_ : geo_divs_;

  // Draw the shapes
  for (const auto &gd : geo_divs) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const Polygon ext_ring = pwh.outer_boundary();
      cairo_move_to(cr, ext_ring[0].x() + 75, ly_ - ext_ring[0].y() + 75);

      // Plot each point in exterior ring
      for (unsigned int i = 1; i < ext_ring.size(); ++i) {
        cairo_line_to(cr, ext_ring[i].x() + 75, ly_ - ext_ring[i].y() + 75);
      }

      // Close the exterior ring
      cairo_close_path(cr);

      // Draw holes
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        cairo_move_to(cr, (*h)[0].x() + 75, ly_ - (*h)[0].y() + 75);
        const size_t hsize = (*h).size();
        for (unsigned int i = 1; i <= hsize; ++i) {
          cairo_line_to(cr, (*h)[i % hsize].x() + 75, ly_ - (*h)[i % hsize].y() + 75);
        }
      }
      if (colors || fill_polygons) {
        if (is_input_target_area_missing(gd.id())) {

          // Fill path with dark gray
          cairo_set_source_rgb(cr, 0.9375, 0.9375, 0.9375);
        } else if (colors) {

          // Get color
          const Color col = color_at(gd.id());

          // Fill path
          cairo_set_source_rgb(cr, col.r, col.g, col.b);
        } else if (fill_polygons) {

          // Fill path with default color
          cairo_set_source_rgb(cr, 0.96, 0.92, 0.70);
        }
        cairo_fill_preserve(cr);
      }

      cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
      cairo_stroke(cr);
    }
  }
}

// Outputs a PNG file
void InsetState::write_cairo_polygons_to_png(
  const std::string &fname,
  const bool fill_polygons,
  const bool colors,
  const bool plot_grid)
{
  const auto filename = fname.c_str();
  cairo_surface_t *surface = cairo_image_surface_create(
    CAIRO_FORMAT_ARGB32,
    static_cast<int>(lx_),
    static_cast<int>(ly_));
  cairo_t *cr = cairo_create(surface);
  write_polygons_to_cairo_surface(cr, fill_polygons, colors, plot_grid);
  cairo_surface_write_to_png(surface, filename);
  cairo_destroy(cr);
  cairo_surface_destroy(surface);
}

// Outputs a PS file
void InsetState::write_cairo_polygons_to_ps(
  const std::string &fname,
  const bool fill_polygons,
  const bool colors,
  const bool plot_grid)
{
  const auto filename = fname.c_str();
  cairo_surface_t *surface = cairo_ps_surface_create(filename, lx_, ly_);
  cairo_t *cr = cairo_create(surface);

  // Add comments
  const std::string title = "%%Title: " + fname;
  cairo_ps_surface_dsc_comment(surface, title.c_str());
  cairo_ps_surface_dsc_comment(
    surface,
    "%%Creator: Michael T. Gastner et al.");
  cairo_ps_surface_dsc_comment(surface, "%%For: Humanity");
  cairo_ps_surface_dsc_comment(surface, "%%Copyright: License CC BY");
  cairo_ps_surface_dsc_comment(surface, "%%Magnification: 1.0000");
  write_polygons_to_cairo_surface(cr, fill_polygons, colors, plot_grid);
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

// TODO: DO WE NEED THIS FUNCTION? WOULD IT NOT MAKE MORE SENSE TO ONLY PRINT
// FILE TYPES INDICATED BY COMMAND-LINE FLAGS?
// Outputs both png and ps files
void InsetState::write_cairo_map(
  const std::string &file_name,
  const bool plot_grid)
{
  const auto png_name = file_name + ".png";
  const auto ps_name = file_name + ".ps";

  // Check whether the has all GeoDivs colored
  const bool has_colors = (colors_size() == n_geo_divs());
  write_cairo_polygons_to_png(png_name, true, has_colors, plot_grid);
  write_cairo_polygons_to_ps(ps_name, true, has_colors, plot_grid);
}


void InsetState::write_grid_heatmap_data(const std::string filename)
{
  unsigned int cell_width = 1;
  const unsigned int resolution = 1;

  InsetState is_copy = (*this);
  is_copy.set_geo_divs(geo_divs_original_);
  auto intersections_with_rays = is_copy.intersec_with_parallel_to('x', resolution);
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
        const double area = grid_cell_area(i, j, cell_width);
        f_csv << i << ", " << j << ", " << area << "\n";
      }
    }
  }
  std::cout << "Grid heatmap data written" << std::endl;
}


// Outputs a SVG/PS file of grid heatmap
void InsetState::write_grid_heatmap_image(
  const std::string filename,
  const bool plot_equal_area_map,
  const bool image_format_ps,
  const bool crop_polygons)
{
  // Whether to draw bar on the cairo surface
  const bool draw_bar = true;

  // Create a cairo surface
  cairo_surface_t *surface;
  image_format_ps
    ? surface = cairo_ps_surface_create(filename.c_str(), lx_ + 150, ly_ + 150)
    : surface = cairo_svg_surface_create(filename.c_str(), lx_, ly_);
  cairo_t *cr = cairo_create(surface);

  // Write header
  if (image_format_ps) {
    write_ps_header(filename, surface);
  }

  // White background
  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_rectangle(cr, 0, 0, lx_ + 150, ly_ + 150);
  cairo_fill(cr);

  // Get inset areas
  const double total_ta = total_target_area();
  const double total_ia = total_inset_area();
  
  // std::cout << "Total target area: " << total_ta << std::endl;
  // std::cout << "Total inset area: " << total_ia << std::endl;
  

  // const Bbox bbox_bar = get_bbox_bar(15, 150);

  // Get the max and min grid cell area points
  Point max_area_cell_point, min_area_cell_point;

  std::tie(max_area_cell_point, min_area_cell_point) =
    max_and_min_grid_cell_area_index(1); // width 1 for cells

  const double max_area_cell_point_area =
    grid_cell_area(max_area_cell_point.x(), max_area_cell_point.y(), 1);
  const double min_area_cell_point_area =
    grid_cell_area(min_area_cell_point.x(), min_area_cell_point.y(), 1);

  // Get the max and min grid cell target area per km
  double max_target_area_per_km = grid_cell_target_area_per_km(
    max_area_cell_point.x(),
    max_area_cell_point.y(),
    total_ta,
    total_ia);
  double min_target_area_per_km = grid_cell_target_area_per_km(
    min_area_cell_point.x(),
    min_area_cell_point.y(),
    total_ta,
    total_ia);

  std::cout << "Max target area per km: " << max_target_area_per_km
                << std::endl;
      std::cout << "Min target area per km: " << min_target_area_per_km
                << std::endl
                << std::endl;
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
  
  std::cout << "Nice numbers: " << std::endl;
  for (auto n : nice_numbers) {
    std::cout << n << ", ";
  }
  std::cout << std::endl;

  std::vector<std::pair<double, double>> major_ticks, minor_ticks;

  std::tie(major_ticks, minor_ticks) = get_ticks(
    5,
    min_target_area_per_km,
    max_target_area_per_km,
    min_area_cell_point_area,
    max_area_cell_point_area,
    nice_numbers);

  // Draw colors
  write_grid_colors_to_cairo_surface(cr, plot_equal_area_map, crop_polygons);

  // Draw polygons without color
  write_polygons_to_cairo_surface(cr, false, false, plot_equal_area_map);

  // trim_grid_heatmap(cr, 20);
  
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
  
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_set_line_width(cr, 1.1);
  
  std::vector<double> mj_ticks;
  std::vector<double> mr_ticks;
  double ratio = 512.0/29.0;
  mj_ticks = {(ratio * 5), (ratio * (5+7)), (ratio * (5+ 7 + 7)), ratio * (5 + 7 + 7 + 7)};
  std::vector<std::string> dates = {"11 Oct", "18 Oct", "25 Oct", "1 Nov"};
  font_size = 10;
  double cur = 0.0;
  for(int i = 1; i < 29; i++) {
    cur += ratio;
    mr_ticks.push_back(cur);
  }
  
  // Draw major ticks
  for(int i = 0; i < mj_ticks.size(); i++) {
    double tick = mj_ticks[i];
    auto date = dates[i];
    cairo_move_to(cr, tick + 75, ly_ - half_tick_width + 75);
    cairo_line_to(cr, tick + 75, ly_ + half_tick_width + 75);
    cairo_fill_preserve(cr);
    cairo_stroke(cr);
    
    cairo_text_extents_t extents;
    cairo_text_extents(cr, date.c_str(), &extents);
    double move = extents.width;
    cairo_move_to(cr, tick + 90 - move, ly_ + half_tick_width * 2 + 85);
    cairo_show_text(cr, date.c_str());
  }
  
  cairo_set_line_width(cr, 0.8);
  // Draw minor ticks
  for(int i = 0; i < mr_ticks.size(); i++) {
    double tick = mr_ticks[i];
    cairo_move_to(cr, tick + 75, ly_ - half_tick_width/2 + 75);
    cairo_line_to(cr, tick + 75, ly_ + half_tick_width/2 + 75);
    cairo_fill_preserve(cr);
    cairo_stroke(cr);
  }
  
  std::vector<double> mj_ticks_y_axis;
  std::vector<double> mr_ticks_y_axis;
  ratio = 512.0/10.0;
  mj_ticks_y_axis = {0, ratio * 2, ratio * 4, ratio * 6, ratio * 8, ratio * 10};
  std::vector<double> dates_y_axis = {0, 20, 40, 60, 80, 100};
  mr_ticks_y_axis = {ratio * 1, ratio * 2, ratio * 3, ratio * 4, ratio * 5, ratio * 6, ratio * 7, ratio * 8, ratio * 9, ratio * 10};
  
  cairo_set_line_width(cr, 1.1);
  
  // Draw major ticks
  for(int i = 0; i < mj_ticks_y_axis.size(); i++) {
    double tick = mj_ticks_y_axis[i];
    auto date = dates_y_axis[i];
    cairo_move_to(cr, 75 - half_tick_width, ly_ + 75 - tick);
    cairo_line_to(cr, 75 + half_tick_width, ly_ + 75 - tick);
    cairo_fill_preserve(cr);
    cairo_stroke(cr);
    
    right_aligned_text(cr, date, 75 - half_tick_width, ly_ + 73 - tick + font_size/2, font_size);
  }
  
  cairo_set_line_width(cr, 0.8);
  // Draw minor ticks
  for(int i = 0; i < mr_ticks_y_axis.size(); i++) {
    double tick = mr_ticks_y_axis[i];
    cairo_move_to(cr, 75 - half_tick_width/2, ly_ + 75 - tick);
    cairo_line_to(cr, 75 + half_tick_width/2, ly_ + 75 - tick);
    cairo_fill_preserve(cr);
    cairo_stroke(cr);
  }

  if (draw_bar) {
    std::string bar_filename = "bar_" + filename;

    // Create a cairo bar_surface
    cairo_surface_t *bar_surface =
      (image_format_ps
         ? cairo_ps_surface_create(bar_filename.c_str(), 160, 400)
         : cairo_svg_surface_create(bar_filename.c_str(), 160, 400));

    cairo_t *bar_cr = cairo_create(bar_surface);

    // Write header
    if (image_format_ps) {
      write_ps_header(bar_filename, bar_surface);
    }
  
    std::cout << min_area_cell_point_area << " " << max_area_cell_point_area
              << std::endl;
              
    // Write bar
    write_grid_heatmap_bar_to_cairo_surface(
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
void write_density_bar_to_cairo_surface(
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
    std::max(max_value, 1.0);

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

// Adds a square legend outside the map
void InsetState::write_legend_to_cairo_surface(
  cairo_t *cr,
  bool equal_area_map)
{
  std::pair<double, unsigned int> legend_info =
    equal_area_map ? get_km_legend_length()
                   : get_visual_variable_legend_length();

  Bbox bb = bbox();
  Point legend_pos(bb.xmin() * 0.5, bb.ymin() * 0.5);
  double legend_length = legend_info.first;
  unsigned int legend_value = legend_info.second;
  cairo_set_line_width(cr, 1);
  Color color("black");
  cairo_set_source_rgb(cr, color.r, color.g, color.b);
  cairo_rectangle(
    cr,
    legend_pos.x(),
    legend_pos.y(),
    legend_length,
    legend_length);
  cairo_stroke(cr);

  double x = legend_pos.x() + (legend_length * 1.25);
  double y = legend_pos.y() + (legend_length * 0.5);

  unsigned int font_size = 12;
  cairo_select_font_face(
    cr,
    "Sans",
    CAIRO_FONT_SLANT_NORMAL,
    CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size(cr, font_size);
  std::string legend_label = std::to_string(legend_value);

  // Add known units
  if (equal_area_map) {
    legend_label += " km²";
  }

  cairo_move_to(cr, x, y + (font_size / 2.0));
  cairo_show_text(cr, legend_label.c_str());
}

void InsetState::write_density_image(
  const std::string filename,
  const double *density,
  const bool plot_pycnophylactic,
  const bool image_format_ps)
{
  // Whether to draw bar on the cairo surface
  const bool draw_bar = true;
  cairo_surface_t *surface =
    (image_format_ps ? cairo_ps_surface_create(filename.c_str(), lx_ + 150, ly_ + 150)
                     : cairo_svg_surface_create(filename.c_str(), lx_, ly_));

  // Create a cairo surface
  cairo_t *cr = cairo_create(surface);

  // Write header
  if (image_format_ps) {
    write_ps_header(filename, surface);
  }

  // White background
  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_rectangle(cr, 0, 0, lx_ + 150, ly_ + 150);
  cairo_fill(cr);

  cairo_set_line_width(cr, 0);

  // Determine range of densities
  const double dens_min = dens_min_;
  const double dens_mean = dens_mean_;
  const double dens_max = dens_max_;
  const double each_grid_cell_area_km = grid_cell_area_km(0, 0);
  std::cout << std::endl << "Each grid cell area: " << each_grid_cell_area_km << " km²" << std::endl;
  std::cout << "Minimum density: " << dens_min << " km²" << std::endl;
  std::cout << "Mean density: " << dens_mean << " km²" << std::endl;
  std::cout << "Maximum density: " << dens_max << " km²" << std::endl;
  
  // Crop it too
  if (plot_pycnophylactic) {
    unsigned int cell_width = 1;

    // Clip to shape and print in sequential scale
    for (const auto &gd : geo_divs_) {
      for (const auto &pwh : gd.polygons_with_holes()) {
        const auto ext_ring = pwh.outer_boundary();
        cairo_new_path(cr);
        cairo_move_to(cr, ext_ring[0].x() + 75, ly_ - ext_ring[0].y() + 75);

        // Plot each point in exterior ring
        for (unsigned int i = 1; i < ext_ring.size(); ++i) {
          cairo_line_to(cr, ext_ring[i].x() + 75, ly_ - ext_ring[i].y() + 75);
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
            Color color = grid_cell_color(target_area_km, 69.9668, 5.15995);

            // Get four points of the square
            double x_min = i - 0.5 * sq_overlap;
            double y_min = j - 0.5 * sq_overlap;
            double x_max = i + 1 + 0.5 * sq_overlap;
            double y_max = j + 1 + 0.5 * sq_overlap;

            cairo_move_to(cr, x_min + 75, ly_ - y_min + 75);
            cairo_line_to(cr, x_max + 75, ly_ - y_min + 75);
            cairo_line_to(cr, x_max + 75, ly_ - y_max + 75);
            cairo_line_to(cr, x_min + 75, ly_ - y_max + 75);

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
    for (unsigned int i = 0; i < lx_; ++i) {
      for (unsigned int j = 0; j < ly_; ++j) {
        Color color =
          heatmap_color(density[i * ly_ + j], dens_min, dens_mean, dens_max);

        // Get four points of the square
        double x_min = i - 0.5 * sq_overlap;
        double y_min = j - 0.5 * sq_overlap;
        double x_max = i + 1 + 0.5 * sq_overlap;
        double y_max = j + 1 + 0.5 * sq_overlap;

        cairo_move_to(cr, x_min + 75, ly_ - y_min + 75);
        cairo_line_to(cr, x_max + 75, ly_ - y_min + 75);
        cairo_line_to(cr, x_max + 75, ly_ - y_max + 75);
        cairo_line_to(cr, x_min + 75, ly_ - y_max + 75);

        cairo_set_source_rgb(cr, color.r, color.g, color.b);
        cairo_fill(cr);
        cairo_set_source_rgb(cr, 0, 0, 0);
        cairo_stroke(cr);
      }
    }
  }
  write_polygons_to_cairo_surface(cr, false, false, false);
  
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
  
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_set_line_width(cr, 1.1);
  
  std::vector<double> major_ticks;
  std::vector<double> minor_ticks;
  double ratio = 512.0/29.0;
  major_ticks = {(ratio * 5), (ratio * (5+7)), (ratio * (5+ 7 + 7)), ratio * (5 + 7 + 7 + 7)};
  std::vector<std::string> dates = {"11 Oct", "18 Oct", "25 Oct", "1 Nov"};
  font_size = 10;
  double cur = 0.0;
  for(int i = 1; i < 29; i++) {
    cur += ratio;
    minor_ticks.push_back(cur);
  }
  
  // Draw major ticks
  for(int i = 0; i < major_ticks.size(); i++) {
    double tick = major_ticks[i];
    auto date = dates[i];
    cairo_move_to(cr, tick + 75, ly_ - half_tick_width + 75);
    cairo_line_to(cr, tick + 75, ly_ + half_tick_width + 75);
    cairo_fill_preserve(cr);
    cairo_stroke(cr);
    
    cairo_text_extents_t extents;
    cairo_text_extents(cr, date.c_str(), &extents);
    double move = extents.width;
    cairo_move_to(cr, tick + 90 - move, ly_ + half_tick_width * 2 + 85);
    cairo_show_text(cr, date.c_str());
  }
  
  cairo_set_line_width(cr, 0.8);
  // Draw minor ticks
  for(int i = 0; i < minor_ticks.size(); i++) {
    double tick = minor_ticks[i];
    cairo_move_to(cr, tick + 75, ly_ - half_tick_width/2 + 75);
    cairo_line_to(cr, tick + 75, ly_ + half_tick_width/2 + 75);
    cairo_fill_preserve(cr);
    cairo_stroke(cr);
  }
  
  std::vector<double> major_ticks_y_axis;
  std::vector<double> minor_ticks_y_axis;
  ratio = 512.0/10.0;
  major_ticks_y_axis = {0, ratio * 2, ratio * 4, ratio * 6, ratio * 8, ratio * 10};
  std::vector<double> dates_y_axis = {0, 20, 40, 60, 80, 100};
  minor_ticks_y_axis = {ratio * 1, ratio * 2, ratio * 3, ratio * 4, ratio * 5, ratio * 6, ratio * 7, ratio * 8, ratio * 9, ratio * 10};
  
  cairo_set_line_width(cr, 1.1);
  
  // Draw major ticks
  for(int i = 0; i < major_ticks_y_axis.size(); i++) {
    double tick = major_ticks_y_axis[i];
    auto date = dates_y_axis[i];
    cairo_move_to(cr, 75 - half_tick_width, ly_ + 75 - tick);
    cairo_line_to(cr, 75 + half_tick_width, ly_ + 75 - tick);
    cairo_fill_preserve(cr);
    cairo_stroke(cr);
    
    right_aligned_text(cr, date, 75 - half_tick_width, ly_ + 73 - tick + font_size/2, font_size);
  }
  
  cairo_set_line_width(cr, 0.8);
  // Draw minor ticks
  for(int i = 0; i < minor_ticks_y_axis.size(); i++) {
    double tick = minor_ticks_y_axis[i];
    cairo_move_to(cr, 75 - half_tick_width/2, ly_ + 75 - tick);
    cairo_line_to(cr, 75 + half_tick_width/2, ly_ + 75 - tick);
    cairo_fill_preserve(cr);
    cairo_stroke(cr);
  }
  
  
  if (draw_bar && !plot_pycnophylactic) {
    std::string bar_filename = "bar_" + filename;

    // Create a cairo bar_surface
    cairo_surface_t *bar_surface =
      (image_format_ps
         ? cairo_ps_surface_create(bar_filename.c_str(), 160, 400)
         : cairo_svg_surface_create(bar_filename.c_str(), 160, 400));

    cairo_t *bar_cr = cairo_create(bar_surface);

    // Write header
    if (image_format_ps) {
      write_ps_header(bar_filename, bar_surface);
    }

    write_density_bar_to_cairo_surface(
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

void InsetState::write_intersections_image(
  unsigned int res,
  const bool image_format_ps)
{
  std::string filename = inset_name() + "_intersections_" +
                         std::to_string(n_finished_integrations());

  // Update extension
  image_format_ps ? filename += ".ps" : filename += ".svg";

  // Calculating intersections
  std::vector<Segment> intersections = intersecting_segments(res);
  cairo_surface_t *surface;

  // Create a cairo surface
  image_format_ps
    ? surface = cairo_ps_surface_create(filename.c_str(), lx_, ly_)
    : surface = cairo_svg_surface_create(filename.c_str(), lx_, ly_);
  cairo_t *cr = cairo_create(surface);

  // Write header
  if (image_format_ps) {
    write_ps_header(filename, surface);
  }

  write_polygons_to_cairo_surface(cr, false, false, false);

  cairo_set_line_width(cr, 0.0001 * std::min(lx_, ly_));

  for (auto seg : intersections) {
    // Move to starting coordinates
    cairo_move_to(cr, seg[0][0], ly_ - seg[0][1]);

    // Draw line
    cairo_line_to(cr, seg[1][0], ly_ - seg[1][1]);

    // line with red and stroke
    cairo_set_source_rgb(cr, 1, 0, 0);
    cairo_stroke(cr);
  }
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

double get_font_size(
  cairo_t *cr,
  const char *label,
  const Point label_pt,
  const GeoDiv gd)
{
  for (double fsize = max_font_size; fsize >= min_font_size; fsize -= 0.5) {
    cairo_set_font_size(cr, fsize);
    cairo_text_extents_t extents;
    cairo_text_extents(cr, label, &extents);
    const Polygon_with_holes largest_pwh = gd.largest_polygon_with_holes();

    // Bounding box of the label
    const CGAL::Bbox_2 bb(
      label_pt.x() - 0.5 * extents.width,
      label_pt.y() - 0.5 * extents.height,
      label_pt.x() + 0.5 * extents.width,
      label_pt.y() + 0.5 * extents.height);

    // Vector of bounding-box edge points
    std::vector<Point> bb_edge_points;
    for (unsigned int i = 0; i <= 1; ++i) {
      for (unsigned int j = 0; j <= 5; ++j) {
        bb_edge_points.push_back(Point(
          (j * bb.xmin() + (5 - j) * bb.xmax()) / 5,
          (i * bb.ymin() + (1 - i) * bb.ymax())));
      }
    }
    if (all_points_inside_exterior_ring(bb_edge_points, largest_pwh)) {
      return fsize;
    }
  }
  return 0.0;
}

void InsetState::write_labels_to_cairo_surface(cairo_t *cr)
{
  for (const auto &gd : geo_divs()) {
    const std::string label = label_at(gd.id());
    const char *const label_char = label.c_str();
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_select_font_face(
      cr,
      "sans-serif",
      CAIRO_FONT_SLANT_NORMAL,
      CAIRO_FONT_WEIGHT_NORMAL);
    const Point label_pt = gd.point_on_surface_of_geodiv();

    // Get size of label
    const double fsize = get_font_size(cr, label_char, label_pt, gd);
    cairo_text_extents_t extents;

    // Draw label only if appropriate size is found
    if (fsize > 0.0) {
      cairo_set_font_size(cr, fsize);
      cairo_text_extents(cr, label_char, &extents);
      const double x = label_pt.x() - (extents.width / 2 + extents.x_bearing);
      const double y =
        ly_ - label_pt.y() - (extents.height / 2 + extents.y_bearing);
      cairo_move_to(cr, x, y);
      cairo_show_text(cr, label_char);
    }
  }
}

// Returns edge points of a grid cell and treats
// them as a polygon.
Polygon InsetState::grid_cell_edge_points(
  unsigned int x,
  unsigned int y,
  unsigned int cell_width = plotted_cell_length,
  bool plot_equal_area_map = false)
{
  Polygon cell_edge_points;
  const boost::multi_array<XYPoint, 2> &proj =
    plot_equal_area_map ? identity_proj_ : cum_proj_;

  // Horizontal lower edge points
  for (unsigned int i = x; i < x + cell_width; ++i) {
    double x_coor_trans = proj[i][y].x;
    double y_coor_trans = proj[i][y].y;
    cell_edge_points.push_back(Point(x_coor_trans, y_coor_trans));
  }

  // Vertical right edge points
  for (unsigned int i = y; i < y + cell_width; ++i) {
    double x_coor_trans = proj[x + cell_width][i].x;
    double y_coor_trans = proj[x + cell_width][i].y;
    cell_edge_points.push_back(Point(x_coor_trans, y_coor_trans));
  }

  // Horizontal upper edge points
  for (unsigned int i = x + cell_width; i > x; --i) {
    double x_coor_trans = proj[i][y + cell_width].x;
    double y_coor_trans = proj[i][y + cell_width].y;
    cell_edge_points.push_back(Point(x_coor_trans, y_coor_trans));
  }

  // Vertical left edge points
  for (unsigned int i = y + cell_width; i > y; --i) {
    double x_coor_trans = proj[x][i].x;
    double y_coor_trans = proj[x][i].y;
    cell_edge_points.push_back(Point(x_coor_trans, y_coor_trans));
  }

  // Complete the polygon by making first and last point same
  cell_edge_points.push_back(cell_edge_points[0]);

  return cell_edge_points;
}

// Returns grid cell area based on edge points
double InsetState::grid_cell_area(
  unsigned int x,
  unsigned int y,
  unsigned int cell_width)
{
  // Taking absolule to ensure we get the area irrespective of direction
  return abs(grid_cell_edge_points(x, y, cell_width).area());
}

// Returns the largest and smallest grid cell area to be used for
// grid heatmap generation
std::pair<double, double> InsetState::max_and_min_grid_cell_area(
  unsigned int cell_width)
{
  // Initialize max and min area
  double max_area = -dbl_inf;
  double min_area = dbl_inf;

  // Iterate over grid cells
  for (unsigned int i = 0; i < lx_ - cell_width; i += cell_width) {
    for (unsigned int j = 0; j < ly_ - cell_width; j += cell_width) {
      const double area = grid_cell_area(i, j, cell_width);
      max_area = std::max(max_area, area);
      min_area = std::min(min_area, area);
    }
  }
  return std::make_pair(max_area, min_area);
}

// Returns the index of largest and smallest grid cell area to be used for
// grid heatmap generation
std::pair<Point, Point> InsetState::max_and_min_grid_cell_area_index(
  unsigned int cell_width)
{
  // Initialize max and min area
  double max_area = -dbl_inf;
  double min_area = dbl_inf;
  unsigned int max_i = 0;
  unsigned int max_j = 0;
  unsigned int min_i = 0;
  unsigned int min_j = 0;

  // Iterate over grid cells
  for (unsigned int i = 0; i < lx_ - cell_width; i += cell_width) {
    for (unsigned int j = 0; j < ly_ - cell_width; j += cell_width) {
      const double area = grid_cell_area(i, j, cell_width);
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

    // // Turqoise to brown
    Color("#543005"),
    Color("#8c510a"),
    Color("#bf812d"),
    Color("#dfc27d"),
    Color("#f6e8c3"),
    Color("#ffffff"),
    Color("#c7eae5"),
    Color("#80cdc1"),
    Color("#35978f"),
    Color("#01665e"),
    Color("#003c30")

    // // Original
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
    // Color("#ffffff"),
    // Color("#efedf5"),
    // Color("#dadaeb"),
    // Color("#bcbddc"),
    // Color("#9e9ac8"),
    // Color("#807dba"),
    // Color("#6a51a3"),
    // Color("#54278f"),
    // Color("#3f007d")

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
    Color(1.000, 0.961, 0.941),
    Color(0.996, 0.878, 0.824),
    Color(0.988, 0.733, 0.631),
    Color(0.988, 0.572, 0.447),
    Color(0.984, 0.416, 0.290),
    Color(0.937, 0.231, 0.173),
    Color(0.796, 0.094, 0.114),
    Color(0.647, 0.058, 0.082),
    Color(0.404, 0.000, 0.050)
  };
  int n_categories = colors.size();

  // Normalize area to [0,1] and make it logarithmic
  double ratio = ((area) - (min_area)) / ((max_area) - (min_area));

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
std::vector<std::vector<Color>> InsetState::grid_cell_colors(
  unsigned int cell_width)
{

  // Initialize max and min area
  double max_area, min_area;
  std::tie(max_area, min_area) = max_and_min_grid_cell_area(cell_width);

  // Initialize colors
  std::vector<std::vector<Color>> colors(
    lx_ - cell_width,
    std::vector<Color>(ly_ - cell_width));

  // Iterate over grid cells
  for (unsigned int i = 0; i < lx_ - cell_width; i += cell_width) {
    for (unsigned int j = 0; j < ly_ - cell_width; j += cell_width) {
      const double area = grid_cell_area(i, j, cell_width);
      colors[i][j] = grid_cell_color(area, max_area, min_area);
    }
  }
  return colors;
}

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

void write_grid_heatmap_bar_to_cairo_surface(
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
      (value_at_gradient_segment),
      (max_value),
      (min_value));
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

  std::cout << "Min_value: " << min_value << std::endl;
  // Draw the ticks and nice_numbers
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
  cairo_set_line_width(cr, 1.1);
  for (auto tick : major_ticks) {
    double area = tick.first;
    int NiceNumber = tick.second;
    if (area > min_value and area < max_value) {
      double y =
        (((area) - (min_value)) / ((max_value) - (min_value))) *
          (ymax_bar - ymin_bar) + ymin_bar;
      std::cout << "y: " << y << std::endl;
      std::cout << "NiceNumber: " << NiceNumber << std::endl;
      std::cout << "area: " << area << std::endl;
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
        (((area) - (min_value)) / ((max_value) - (min_value))) *
          (ymax_bar - ymin_bar) +
        ymin_bar;

      cairo_move_to(cr, xmin_bar + half_tick_width, ly - y);
      cairo_line_to(cr, xmin_bar - half_tick_width, ly - y);
      cairo_stroke(cr);
    }
  }
  
    double x = xmin_bar - 30;
    double y = ly - ymax_bar - 24;

    font_size = 10;
    cairo_select_font_face(
      cr,
      "Sans",
      CAIRO_FONT_SLANT_NORMAL,
      CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, font_size);
    std::string legend_label = "Cases per day and";
    std::string legend_label2 = "per year of age";

    cairo_move_to(cr, x, y + (font_size / 2.0));
    cairo_show_text(cr, legend_label.c_str());
    
    cairo_move_to(cr, x + 6.5, y + 2 + (font_size * 1.5));
    cairo_show_text(cr, legend_label2.c_str());
  
  // half_tick_width = bar_width / 8.0;

  // // Write max value at top of bar
  // right_aligned_text(
  //   cr,
  //   std::floor(max_value),
  //   xmin_bar - half_tick_width,
  //   ly - ymax_bar,
  //   font_size);

  // Write mean value
  // right_aligned_text(
  //   cr,
  //   mean_value,
  //   xmin_bar,
  //   ly - ymean_bar + (font_size / 4.0),
  //   font_size);

  // // Write min_value
  // right_aligned_text(
  //   cr,
  //   std::ceil(min_value),
  //   xmin_bar - half_tick_width,
  //   ly - ymin_bar + (font_size * 1.2),
  //   font_size);
}

// Given padding, makes padding area white
void InsetState::trim_grid_heatmap(cairo_t *cr, double padding)
{
  // Canvas dimension
  Bbox is_bbox = bbox();
  double xmin = is_bbox.xmin() - padding;
  double xmax = is_bbox.xmax() + padding;
  double ymin = is_bbox.ymin() - padding;
  double ymax = is_bbox.ymax() + padding;

  // Color white outside is_bbox
  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_rectangle(cr, 0, 0, xmin, ymax);
  cairo_rectangle(cr, xmax, 0, lx_ - xmax, ymax);
  cairo_rectangle(cr, 0, ymax, lx_, ly_ - ymax);
  cairo_rectangle(cr, 0, 0, lx_, ymin);
  cairo_fill(cr);
}

// Writes grid cells and colors them if required
void InsetState::write_grids_to_cairo_surface(cairo_t *cr)
{

  // Set line width of grid lines
  cairo_set_line_width(cr, 5e-4 * std::min(lx_, ly_));
  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);

  // Iterate over grid cells
  for (unsigned int i = 0; i < lx_ - plotted_cell_length;
       i += plotted_cell_length) {
    for (unsigned int j = 0; j < ly_ - plotted_cell_length;
         j += plotted_cell_length) {

      // Draw grid cell by connecting edge points
      const Polygon cell_edge_points = grid_cell_edge_points(i, j);
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
      cairo_stroke(cr);
    }
  }
}

void InsetState::write_grid_colors_to_cairo_surface(
  cairo_t *cr,
  bool plot_equal_area_map,
  bool crop_polygons)
{
  unsigned int cell_width = 1;

  // Get colors
  const auto colors = grid_cell_colors(cell_width);

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
          cell_edge_points[0].x() + 75,
          ly_ - cell_edge_points[0].y() + 75);
        for (unsigned int k = 1; k < cell_edge_points.size(); ++k) {
          cairo_line_to(
            cr,
            cell_edge_points[k].x() + 75,
            ly_ - cell_edge_points[k].y() + 75);
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
  return;
}

// Outputs a SVG/PS file with polygons, labels, and grids (if required)
void InsetState::write_map_image(
  const std::string filename,
  const bool fill_polygons,
  const bool plot_grid,
  const bool plot_labels,
  const bool image_format_ps,
  const bool equal_area_map)
{

  // Check whether the map has all GeoDivs colored
  const bool colors = (colors_size() == n_geo_divs());
  cairo_surface_t *surface;

  // Create a cairo surface
  image_format_ps
    ? surface = cairo_ps_surface_create(filename.c_str(), lx_, ly_)
    : surface = cairo_svg_surface_create(filename.c_str(), lx_, ly_);
  cairo_t *cr = cairo_create(surface);

  // Write header
  if (image_format_ps) {
    write_ps_header(filename, surface);
  }

  // Draw polygons with color
  write_polygons_to_cairo_surface(cr, fill_polygons, colors, false);

  // Place labels
  if (plot_labels) {
    write_labels_to_cairo_surface(cr);
  }

  // Add legend
  write_legend_to_cairo_surface(cr, equal_area_map);

  // Draw grid without color
  if (plot_grid) {
    write_grids_to_cairo_surface(cr);
  }
  cairo_show_page(cr);
  cairo_surface_destroy(surface);
  cairo_destroy(cr);
}

// Given coordinates in 512 by 512 coordinate system, returns the corresponding
// coordinates in the albers projection coordinate system
Polygon InsetState::transform_to_albers_coor(Polygon edge_points)
{
  Transformation scale(CGAL::SCALING, latt_const_);

  Polygon cell_edge_points_albers = transform(scale, edge_points);
  return cell_edge_points_albers;
}

// Given area in the albers projection coordinate system, returns the
// corresponding area in the square km^2
double albers_area_to_earth_area(const double albers_area)
{
  return (albers_area * earth_surface_area) / (4 * pi);
}

double InsetState::grid_cell_area_km(
  const unsigned int i = 0,
  const unsigned int j = 0)
{
  Polygon cell_edge_points = grid_cell_edge_points(i, j, 1, true);
  // const Polygon cell_edge_points_albers =
  //   transform_to_albers_coor(cell_edge_points);
  // const double cell_area = cell_edge_points_albers.area();
  // const double cell_area_km = albers_area_to_earth_area(cell_area);
  double cell_area = cell_edge_points.area();
  double ratio = 2900.0/ (512.0*512.0);
  std::cout << "ratio: " << ratio << std::endl;
  std::cout << "cell_area: " << cell_area << std::endl;
  double cell_area_per_unit = cell_area * ratio;
  return cell_area_per_unit;
}

std::pair<double, unsigned int> InsetState::get_km_legend_length()
{

  return std::pair<double, unsigned int>(0, 0);

  // 1% of the total area, rounded up to the nearest power of 10
  double min_length = 0.02 * ((lx() + ly()) / 2.0);
  double max_length = 0.06 * ((lx() + ly()) / 2.0);
  double unit_square_area = grid_cell_area_km();
  unsigned int legend_area = total_inset_area() * 0.01;
    // pow(10.0, ceil(std::log10());
  double length = sqrt(legend_area / unit_square_area);

  while (length < min_length) {
    length *= 2;
    legend_area *= 4;
  }

  while (length > max_length) {
    length /= 2;
    legend_area /= 4;
  }

  return std::pair<double, unsigned int>(length, legend_area);
}

double InsetState::grid_cell_target_area(
  const unsigned int i,
  const unsigned int j,
  const double total_target_area,
  const double total_inset_area)
{
  const Polygon cell_edge_points = grid_cell_edge_points(i, j, 1, false);
  const double cell_area = cell_edge_points.area();

  const double cell_target_area =
    (cell_area * total_target_area) / total_inset_area;

  return cell_target_area;
}

std::pair<double, unsigned int> InsetState::get_visual_variable_legend_length()
{
  unsigned int legend_area = total_target_area() * 0.01;
    // pow(10.0, ceil(std::log10());

  double unit_square_area = total_target_area() / total_inset_area();
  double min_length = 0.02 * ((lx() + ly()) / 2.0);
  double max_length = 0.06 * ((lx() + ly()) / 2.0);
  double length = sqrt(legend_area / unit_square_area);

  while (length < min_length) {
    length *= 2;
    legend_area *= 4;
  }

  while (length > max_length) {
    length /= 2;
    legend_area /= 4;
  }

  return std::pair<double, unsigned int>(length, legend_area);
}

double InsetState::grid_cell_target_area_per_km(
  const unsigned int i,
  const unsigned int j,
  const double total_target_area,
  const double total_inset_area)
{
  const double cell_target_area =
    grid_cell_target_area(i, j, total_target_area, total_inset_area);
  const double cell_area_km = grid_cell_area_km(i, j);

  const double cell_target_area_per_km = cell_target_area / cell_area_km;

  return cell_target_area_per_km;
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
    std::cout << "area: " << area << std::endl;
    std::cout << "niceNumber: " << niceNumber << std::endl;
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

Bbox InsetState::get_bbox_bar(const double bar_width, const double bar_height)
{
  const Bbox is_bbox = bbox();

  // Position the bar 25 pixels to the right of the is_bbox
  const double xmin_bar = is_bbox.xmax() + 35;
  const double xmax_bar = xmin_bar + bar_width;

  // Position the bar at the middle of the is_bbox y coordinates
  double ymid_bar = (is_bbox.ymax() + is_bbox.ymin()) / 2 - 25;
  double ymin_bar = ymid_bar - bar_height / 2;
  double ymax_bar = ymid_bar + bar_height / 2;

  return Bbox(xmin_bar, ymin_bar, xmax_bar, ymax_bar);
}

std::vector<int> get_nice_numbers_for_bar(const double max_target_area_per_km)
{
  std::vector<int> nice_numbers;
  int NiceNumber = 0;
  nice_numbers.push_back(NiceNumber);
  while (NiceNumber < max_target_area_per_km) {
    NiceNumber = NiceNumber + 10;
    nice_numbers.push_back(NiceNumber);
  }
  return nice_numbers;
}
