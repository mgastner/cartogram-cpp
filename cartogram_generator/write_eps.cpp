#include "constants.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include <fstream>

void write_eps_header_and_definitions(std::ofstream &eps_file,
                                      std::string eps_name,
                                      InsetState *inset_state)
{
  // Current time
  time_t tt;
  time(&tt);
  struct tm *ti = localtime(&tt);

  // Header
  eps_file << "%!PS-Adobe-3.0 EPSF-3.0\n";
  eps_file << "%%Title: " << eps_name << "\n";
  eps_file << "%%Creator: Michael T. Gastner et al.\n";
  eps_file << "%%For: Humanity\n";
  eps_file << "%%Copyright: License CC BY-NC-ND 2.0\n";
  eps_file << "%%CreationDate: " << asctime(ti);
  eps_file << "%%Pages: 1\n";
  eps_file << "%%BoundingBox: 0 0 "
           << inset_state->lx() << " "
           << inset_state->ly() << "\n";
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

void write_polygons_to_eps(std::ofstream &eps_file,
                           bool fill_polygons,
                           bool colors,
                           bool plot_graticule,
                           InsetState *inset_state)
{
  eps_file << 0.001 * std::min(inset_state->lx(), inset_state->ly())
           << " slw\n";
  for (auto gd : inset_state->geo_divs()) {
    for (auto pwh : gd.polygons_with_holes()) {
      Polygon ext_ring = pwh.outer_boundary();

      // Move to starting coordinates
      eps_file << "n " << ext_ring[0][0] << " " << ext_ring[0][1] << " m\n";

      // Plot each point in exterior ring
      for (unsigned int i = 1; i < ext_ring.size(); ++i) {
        eps_file << ext_ring[i][0] << " " << ext_ring[i][1] << " l\n";
      }

      // Close path
      eps_file << "c\n";

      // Plot holes
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
        Polygon hole = *hci;
        eps_file << hole[0][0] << " " << hole[0][1] << " m\n";
        for (unsigned int i = 1; i < hole.size(); ++i) {
          eps_file << hole[i][0] << " " << hole[i][1] << " l\n";
        }
        eps_file << "c\n";
      }
      if (colors || fill_polygons) {

        // Save path before filling it
        eps_file << "gsave\n";

        // Check if target area was initially was missing
        if (inset_state->is_input_target_area_missing(gd.id())) {

          // Fill path with dark-grey
          eps_file << "0.9375 0.9375 0.9375 srgb f\n";

        } else if (colors) {

          // Get color
          Color col = inset_state->colors_at(gd.id());

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

  if (plot_graticule) {
    boost::multi_array<XYPoint, 2> &cum_proj =
      *inset_state->ref_to_cum_proj();
    unsigned int graticule_line_spacing = 7;

    // Set line width of graticule lines
    eps_file << 0.0005 * std::min(inset_state->lx(), inset_state->ly())
             << " slw\n";

    // Vertical graticule lines
    for (unsigned int i = 0;
         i <= inset_state->lx();
         i += graticule_line_spacing) {
      eps_file << cum_proj[i][0].x << " " << cum_proj[i][0].y << " m\n";
      for (unsigned int j = 1; j < inset_state->ly(); ++j) {
        eps_file << cum_proj[i][j].x << " " << cum_proj[i][j].y << " l\n";
      }
      eps_file << "s\n";
    }

    // Horizontal graticule lines
    for (unsigned int j = 0;
         j <= inset_state->ly();
         j += graticule_line_spacing) {
      eps_file << cum_proj[0][j].x << " " << cum_proj[0][j].y << " m\n";
      for (unsigned int i = 1; i < inset_state->lx(); ++i) {
        eps_file << cum_proj[i][j].x << " " << cum_proj[i][j].y << " l\n";
      }
      eps_file << "s\n";
    }
  }
  return;
}

void write_intersections_to_eps(std::ofstream &eps_file,
                                std::vector<Polygon_with_holes> intersections)
{
  for (auto pwh : intersections) {
    Polygon ext_ring = pwh.outer_boundary();

    // Move to starting coordinates
    eps_file << "n " << ext_ring[0][0] << " " << ext_ring[0][1] << " m\n";

    // Plot each point in exterior ring
    for (unsigned int i = 1; i < ext_ring.size(); ++i) {
      eps_file << ext_ring[i][0] << " " << ext_ring[i][1] << " l\n";
    }

    // Close path
    eps_file << "c\n";

    // Plot holes
    for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
      Polygon hole = *hci;
      eps_file << hole[0][0] << " " << hole[0][1] << " m\n";
      for (unsigned int i = 1; i < hole.size(); ++i) {
        eps_file << hole[i][0] << " " << hole[i][1] << " l\n";
      }
      eps_file << "c\n";
    }

    // Fill path with red
    eps_file << "1 0 0 srgb f\n";
  }
}

void write_map_to_eps(std::string eps_name, bool plot_graticule,
                      InsetState *inset_state)
{
  std::ofstream eps_file(eps_name);
  write_eps_header_and_definitions(eps_file, eps_name, inset_state);

  // Check whether the has all GeoDivs colored
  bool has_colors = (inset_state->colors_size() == inset_state->n_geo_divs());
  write_polygons_to_eps(eps_file,
                        true,
                        has_colors,
                        plot_graticule,
                        inset_state);
  write_intersections_to_eps(eps_file,
                             inset_state->intersections());
  eps_file << "showpage\n";
  eps_file << "%%EOF\n";
  eps_file.close();
  return;
}

// Functions to show a scalar field called "density" as a heat map
double interpolate_for_heatmap(double x,
                               double xmin,
                               double xmax,
                               double ymin,
                               double ymax)
{
  return ((x - xmin) * ymax + (xmax - x) * ymin) / (xmax - xmin);
}

void heatmap_color(double dens,
                   double dens_min,
                   double dens_mean,
                   double dens_max,
                   double *r,
                   double *g,
                   double *b)
{

  // Assign possible categories for red, green, blue
  double red[] = {
    0.33, 0.55, 0.75, 0.87, 0.96, 0.99, 0.78, 0.50, 0.21, 0.00, 0.00
  };
  double green[] = {
    0.19, 0.32, 0.51, 0.76, 0.91, 0.96, 0.92, 0.80, 0.59, 0.40, 0.24
  };
  double blue[] = {
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

void write_density_to_eps(std::string eps_name,
                          double *density,
                          InsetState *inset_state)
{
  std::ofstream eps_file(eps_name);
  write_eps_header_and_definitions(eps_file, eps_name, inset_state);

  // Determine range of densities
  double dens_min = dbl_inf;
  double dens_mean = 0.0;
  double dens_max = -dbl_inf;
  unsigned int n_grid_cells = inset_state->lx() * inset_state->ly();
  for (unsigned int k = 0; k < n_grid_cells; ++k) {
    dens_min = std::min(density[k], dens_min);
    dens_mean += density[k];
    dens_max = std::max(density[k], dens_max);
  }
  dens_mean /= n_grid_cells;
  for (unsigned int i = 0; i < inset_state->lx(); ++i) {
    for (unsigned int j = 0; j < inset_state->ly(); ++j) {
      double r, g, b;
      heatmap_color(density[i*inset_state->ly() + j],
                    dens_min,
                    dens_mean,
                    dens_max,
                    &r, &g, &b);
      eps_file << i - 0.5*sq_overlap << " " << j - 0.5*sq_overlap  << " sq ";
      eps_file << r << " " << g << " " << b << " srgb f\n";
    }
  }
  write_polygons_to_eps(eps_file, false, false, false, inset_state);
  eps_file << "showpage\n";
  eps_file << "%%EOF\n";
  eps_file.close();
  return;
}
