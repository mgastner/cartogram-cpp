// #include "constants.h"
// #include <fstream>

void write_eps_header_and_definitions(std::ofstream &eps_file,
                                      const std::string eps_name,
                                      const unsigned int lx,
                                      const unsigned int ly)
{
  // Current time
  time_t tt;
  time(&tt);
  const struct tm *ti = localtime(&tt);

  // Header
  eps_file << "%!PS-Adobe-3.0 EPSF-3.0\n";
  eps_file << "%%Title: " << eps_name << "\n";
  eps_file << "%%Creator: Michael T. Gastner et al.\n";
  eps_file << "%%For: Humanity\n";
  eps_file << "%%Copyright: License CC BY-NC-ND 2.0\n";
  eps_file << "%%CreationDate: " << asctime(ti);
  eps_file << "%%Pages: 1\n";
  eps_file << "%%BoundingBox: 0 0 "
           << lx << " "
           << ly << "\n";
  eps_file << "%%Magnification: 1.0000\n";
  eps_file << "%%EndComments\n";

  // Definitions
  eps_file << "/m {moveto} def\n";
  eps_file << "/l {lineto} def\n";
  eps_file << "/rl {rlineto} def\n";
  eps_file << "/s {stroke} def\n";
  eps_file << "/n {newpath} def\n";
  eps_file << "/c {closepath} def\n";
  eps_file << "/f {fill} def\n";
  eps_file << "/slw {setlinewidth} def\n";
  eps_file << "/sgry {setgray} def\n";
  eps_file << "/srgb {setrgbcolor} def\n";
  eps_file << "/sq {\n";
  eps_file << "  /ry exch def\n";
  eps_file << "  /rx exch def\n";
  eps_file << "  n rx ry m "
           << 1.0 + sq_overlap
           << " 0 rl 0 "
           << 1.0 + sq_overlap
           << " rl "
           << -1.0 - sq_overlap
           << " 0 rl c\n";
  eps_file << "} def %% square\n";
  return;
}

void InsetState::write_polygons_to_eps(std::ofstream &eps_file,
                                       const bool fill_polygons,
                                       const bool colors)
{
  eps_file << 0.001 * std::min(lx_, ly_)
           << " slw\n";
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const Polygon ext_ring = pwh.outer_boundary();

      // Move to starting coordinates
      eps_file << "n " << ext_ring[0][0] << " " << ext_ring[0][1] << " m\n";

      // Plot each point in exterior ring
      for (unsigned int i = 1; i < ext_ring.size(); ++i) {
        eps_file << ext_ring[i][0] << " " << ext_ring[i][1] << " l\n";
      }

      // Close path
      eps_file << "c\n";

      // Plot holes
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        const Polygon hole = *h;
        eps_file << hole[0][0] << " " << hole[0][1] << " m\n";
        for (unsigned int i = 1; i < hole.size(); ++i) {
          eps_file << hole[i][0] << " " << hole[i][1] << " l\n";
        }
        eps_file << "c\n";
      }
      if (colors || fill_polygons) {

        // Save path before filling it
        eps_file << "gsave\n";

        // Check whether target area was initially missing
        if (is_input_target_area_missing(gd.id())) {

          // Fill path with dark grey
          eps_file << "0.9375 0.9375 0.9375 srgb f\n";

        } else if (colors) {

          // Get color
          Color col = color_at(gd.id());

          // Fill path
          eps_file << col.eps() << "srgb f\n";

        } else if (fill_polygons) {

          // Fill path with default color
          eps_file << "0.96 0.92 0.70 srgb f\n";
        }

        // Restore path
        eps_file << "grestore\n";
      }

      // Stroke path
      eps_file << "0 sgry s\n";
    }
  }
  return;
}

void InsetState::write_graticule_to_eps(std::ofstream &eps_file)
{
  const unsigned int graticule_line_spacing = 7;

  // Set line width of graticule lines
  eps_file << 0.0005 * std::min(lx_, ly_)
           << " slw\n";

  // Vertical graticule lines
  for (unsigned int i = 0;
       i <= lx_;
       i += graticule_line_spacing) {
    eps_file << cum_proj_[i][0].x << " " << cum_proj_[i][0].y << " m\n";
    for (unsigned int j = 1; j < ly_; ++j) {
      eps_file << cum_proj_[i][j].x << " " << cum_proj_[i][j].y << " l\n";
    }
    eps_file << "s\n";
  }

  // Horizontal graticule lines
  for (unsigned int j = 0;
       j <= ly_;
       j += graticule_line_spacing) {
    eps_file << cum_proj_[0][j].x << " " << cum_proj_[0][j].y << " m\n";
    for (unsigned int i = 1; i < lx_; ++i) {
      eps_file << cum_proj_[i][j].x << " " << cum_proj_[i][j].y << " l\n";
    }
    eps_file << "s\n";
  }
  return;
}


void InsetState::write_map_to_eps(const std::string eps_name,
                                  const bool plot_graticule)
{
  std::ofstream eps_file(eps_name);
  write_eps_header_and_definitions(eps_file, eps_name, lx_, ly_);

  // Check whether all GeoDivs are colored
  const bool has_colors =
    (colors_size() == n_geo_divs());
  write_polygons_to_eps(eps_file,
                        true,  // Fill polygons with default color?
                        has_colors);  // Fill polygons with assigned assigned?
  if (plot_graticule) {
    write_graticule_to_eps(eps_file);
  }
  eps_file << "showpage\n";
  eps_file << "%%EOF\n";
  eps_file.close();
  return;
}

// Functions to show a scalar field called "density" as a heat map
double interpolate_for_heatmap(const double x,
                               const double xmin,
                               const double xmax,
                               const double ymin,
                               const double ymax)
{
  return ((x - xmin) * ymax + (xmax - x) * ymin) / (xmax - xmin);
}

void heatmap_color(const double dens,
                   const double dens_min,
                   const double dens_mean,
                   const double dens_max,
                   double *r,
                   double *g,
                   double *b)
{
  // Assign possible categories for red, green, blue
  const double red[] = {
    0.33, 0.55, 0.75, 0.87, 0.96, 0.99, 0.78, 0.50, 0.21, 0.00, 0.00
  };
  const double green[] = {
    0.19, 0.32, 0.51, 0.76, 0.91, 0.96, 0.92, 0.80, 0.59, 0.40, 0.24
  };
  const double blue[] = {
    0.02, 0.04, 0.18, 0.49, 0.76, 0.89, 0.90, 0.76, 0.56, 0.37, 0.19
  };
  double xmin, xmax;
  int color_category;

  // Choose color category
  if (dens > dens_max) {
    *r = red[0];
    *g = green[0];
    *b = blue[0];
    return;
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
    color_category = std::min(color_category, 9);
  } else {
    *r = red[10];
    *g = green[10];
    *b = blue[10];
    return;
  }
  *r = interpolate_for_heatmap(dens,
                               xmin,
                               xmax,
                               red[color_category + 1],
                               red[color_category]);
  *g = interpolate_for_heatmap(dens,
                               xmin,
                               xmax,
                               green[color_category + 1],
                               green[color_category]);
  *b = interpolate_for_heatmap(dens,
                               xmin,
                               xmax,
                               blue[color_category + 1],
                               blue[color_category]);
  return;
}

void InsetState::write_density_to_eps(const std::string eps_name,
                                      const double *density)
{
  std::ofstream eps_file(eps_name);
  write_eps_header_and_definitions(eps_file, eps_name, lx_, ly_);

  // Determine range of densities
  double dens_min = dbl_inf;
  double dens_mean = 0.0;
  double dens_max = -dbl_inf;
  const unsigned int n_grid_cells = lx_ * ly_;
  for (unsigned int k = 0; k < n_grid_cells; ++k) {
    dens_min = std::min(density[k], dens_min);
    dens_mean += density[k];
    dens_max = std::max(density[k], dens_max);
  }
  dens_mean /= n_grid_cells;
  for (unsigned int i = 0; i < lx_; ++i) {
    for (unsigned int j = 0; j < ly_; ++j) {
      double r, g, b;
      heatmap_color(density[i*ly_ + j],
                    dens_min,
                    dens_mean,
                    dens_max,
                    &r, &g, &b);
      eps_file << i - 0.5*sq_overlap << " " << j - 0.5*sq_overlap  << " sq ";
      eps_file << r << " " << g << " " << b << " srgb f\n";
    }
  }
  write_polygons_to_eps(eps_file,
                        false,  // Fill polygons with default color?
                        false);  // Fill polygons with assigned colors?

  eps_file << "showpage\n";
  eps_file << "%%EOF\n";
  eps_file.close();
  return;
}

void InsetState::write_intersections_to_eps(unsigned int res)
{
  std::string eps_name =
    inset_name() +
    "_intersections_" +
    std::to_string(n_finished_integrations()) +
    ".eps";

  // Calculating intersections
  std::vector<Segment> intersections = intersecting_segments(res);

  // Printing intersections to EPS if intersections present
  std::cerr << "Writing " << eps_name << std::endl;
  std::ofstream eps_file(eps_name);
  write_eps_header_and_definitions(eps_file, eps_name, lx_, ly_);
  write_polygons_to_eps(eps_file,
                        false,  // Fill polygons with default color?
                        false);  // Fill polygons with assigned colors?

  // Set line width of intersection lines
  eps_file << 0.0001 * std::min(lx_, ly_)
           << " slw\n";
  for (auto seg : intersections) {

    // Move to starting coordinates
    eps_file << seg[0][0] << " " << seg[0][1] << " m\n";

    // Draw line
    eps_file << seg[1][0] << " " << seg[1][1] << " l\n";

    // Fill line with red and stroke
    eps_file << "1 0 0 srgb s\n";
  }
  eps_file << "showpage\n";
  eps_file << "%%EOF\n";
  eps_file.close();
  return;
}
