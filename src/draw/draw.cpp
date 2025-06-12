#include "common/basic_figures.hpp"
#include "common/canvas.hpp"
#include "common/color.hpp"
#include "common/legend.hpp"
#include "inset_state.hpp"
#include <unordered_set>

void InsetState::write_density_image(const std::string filename)
{
  double *density = rho_init_.as_1d_array();
  std::cerr << "Writing " << filename << std::endl;

  Canvas cvs(filename, lx_, ly_);

  const double dens_min = dens_min_;
  const double dens_mean = dens_mean_;
  const double dens_max = dens_max_;
  const double exterior_density = exterior_density_;

  Color exterior_density_color =
    heatmap_color(exterior_density, dens_min, dens_mean, dens_max);

  cvs.set_fill(exterior_density_color);
  cvs.rectangle(0, 0, lx_, ly_, /*filled=*/true);
  cvs.clear_fill();

  const double stroke_w = 0.0001 * std::min(lx_, ly_);
  cvs.set_stroke(Color{0.0, 0.0, 0.0}, stroke_w);

  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {

      Color color =
        heatmap_color(density[i * ly_ + j], dens_min, dens_mean, dens_max);

      if (color == exterior_density_color)
        continue;

      double x_min = i - 0.5 * sq_overlap;
      double y_min = j - 0.5 * sq_overlap;
      double x_max = i + 1 + 0.5 * sq_overlap;
      double y_max = j + 1 + 0.5 * sq_overlap;

      cvs.set_fill(color);
      cvs.rectangle(
        x_min,
        y_min,
        x_max - x_min,
        y_max - y_min,
        /*filled=*/true);
      cvs.clear_fill();

      cvs.rectangle(
        x_min,
        y_min,
        x_max - x_min,
        y_max - y_min,
        /*filled=*/false);
    }
  }

  write_quadtree_rectangles(cvs, quadtree_bboxes_, Color{0.6, 0.6, 0.6});
  write_polygons(
    cvs,
    /*fill_polygons=*/false,
    /*colours=*/false,
    *this,
    0.025);
}

void InsetState::write_delaunay_triangles(
  const std::string &filename,
  const bool draw_projected_points)
{
  
  
  
   // Skip if redirecting to stdout
  if (args_.redirect_exports_to_stdout) {
    return;
  }
  std::cerr << "Writing " << filename << ".svg" << std::endl;

  Canvas cvs(filename + ".svg", lx_, ly_);

  add_white_background(cvs, *this);

  write_segments(cvs, failed_constraints_dt_projected_, Color{1.0, 0.0, 0.0});
  write_segments(cvs, failed_constraints_, Color{0.0, 0.0, 1.0});

  write_triangles(
    cvs,
    proj_qd_,
    Color{0.6, 0.6, 0.6},
    ly_,
    draw_projected_points);

  write_polygons(
    cvs,
    /*fill_polygons=*/false,
    /*colours       =*/false,
    *this,
    /*line_w =*/0.0,
    Color{1.0, 0.0, 0.0});

  if (draw_projected_points) {
    project_point_set(points_from_densification_);
    project_point_set(points_before_densification_);
  }

  write_point_set(
    cvs,
    points_before_densification_,
    Color{0.0, 0.0, 0.0},
    0.95);
  write_point_set(
    cvs,
    points_from_densification_,
    Color{0.149, 0.545, 0.824},
    0.8);

  failed_constraints_dt_projected_.clear();
  failed_constraints_.clear();
}

void InsetState::write_quadtree(const std::string &filename)
{
   // Skip if redirecting to stdout
  if (args_.redirect_exports_to_stdout) {
    return;
  }
  std::cerr << "Writing " << filename << ".svg" << std::endl;

  Canvas cvs(filename + ".svg", lx_, ly_);

  write_polygons(
    cvs,
    /*fill_polygons=*/false,
    /*colours       =*/false,
    *this);

  write_quadtree_rectangles(cvs, quadtree_bboxes_, Color{0.6, 0.6, 0.6});

  write_polygon_points(cvs, Color{0.0, 0.0, 1.0}, *this);
}

void InsetState::write_intersections_image()
{
  // Skip if redirecting to stdout
  if (args_.redirect_exports_to_stdout) {
    return;
  }
  unsigned int res = intersections_resolution;
  std::string svg_name = inset_name() + "_intersections_" +
                         std::to_string(n_finished_integrations()) + ".svg";

  std::vector<Segment> intersections = intersecting_segments(res);

  Canvas cvs(svg_name, lx_, ly_);

  write_polygons(
    cvs,
    /*fill_polygons=*/false,
    /*colours=*/false,
    *this);

  double w = 0.0001 * std::min(lx_, ly_);
  cvs.set_stroke(Color{1.0, 0.0, 0.0}, w);

  for (const auto &seg : intersections) {
    cvs.move_to(seg[0].x(), seg[0].y());
    cvs.line_to(seg[1].x(), seg[1].y());
    cvs.stroke();
  }
}

void InsetState::write_map(
  const std::string &file_name,
  const bool plot_grid,
  const bool equal_area_map,
  const std::unordered_map<Point, Vector> vectors) const
{
  std::cerr << "Redirect value: " << args_.redirect_exports_to_stdout << std::endl;

  // Skip if redirecting to stdout
  if (args_.redirect_exports_to_stdout) {
    return;
  }
  const auto svg_name = file_name + ".svg";
  std::cerr << "Writing " << svg_name << std::endl;

  Canvas cvs(svg_name, lx(), ly());

  const bool has_colors = (colors_size() == n_geo_divs());
  write_polygons(
    cvs,
    /*fill_polygons=*/true,
    /*colours       =*/has_colors,
    *this);

  if (plot_grid) {
    write_grid(cvs, *this);
    write_legend(cvs, equal_area_map, *this);
  }

  write_vectors(cvs, vectors, lx(), ly());
}
