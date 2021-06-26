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
                           InsetState *inset_state)
{
  eps_file << 0.001 * std::min(inset_state->lx(), inset_state->ly()) << " slw\n";
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
      if (colors) {

        // Getting color
        Color col = inset_state->colors_at(gd.id());

        // Save path before filling it
        eps_file << "gsave\n";

        // Fill path
        eps_file << col.eps() << "srgb f\n";

        // Restore path.
        eps_file << "grestore\n";
      }
      else if (fill_polygons) {

        // Save path before filling it
        eps_file << "gsave\n";

        // Fill path
        eps_file << "0.96 0.92 0.70 srgb f\n";

        // Restore path.
        eps_file << "grestore\n";
      }

      // Stroke path
      eps_file << "0 sgry s\n";
    }
  }
  return;
}

void write_map_to_eps(std::string eps_name, InsetState *inset_state)
{
  std::ofstream eps_file(eps_name);
  write_eps_header_and_definitions(eps_file, eps_name, inset_state);
  write_polygons_to_eps(eps_file,
                        true,
                        // false,
                        !(inset_state->colors_empty()),
                        inset_state);
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
  if (dens > dens_max) {
    *r = red[0];
    *g = green[0];
    *b = blue[0];
    return;
  } else if (dens > dens_mean) {
    color_category = 5 * (dens_max - dens) / (dens_max - dens_mean);
    xmax = dens_max - 0.2 * color_category * (dens_max - dens_mean);
    xmin = xmax - 0.2 * (dens_max - dens_mean);
  } else if (dens > dens_min) {
    color_category = 5 * (dens_mean - dens) / (dens_mean - dens_min) + 5;
    xmax = dens_mean - 0.2 * (color_category - 5) * (dens_mean - dens_min);
    xmin = xmax - 0.2 * (dens_mean - dens_min);
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
                          InsetState *inset_state)
{
  std::ofstream eps_file(eps_name);
  write_eps_header_and_definitions(eps_file, eps_name, inset_state);
  unsigned int lx = inset_state->lx();
  unsigned int ly = inset_state->ly();
  Real2dArray density = inset_state->rho_init();

  // Determine range of density
  double dens_min = density(0, 0);
  double dens_mean = 0.0;
  double dens_max = density(0, 0);
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      dens_min = std::min(density(i, j), dens_min);
      dens_mean += density(i, j);
      dens_max = std::max(density(i, j), dens_max);
    }
  }
  dens_mean /= lx * ly;
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      double r, g, b;
      heatmap_color(density(i, j),
                    dens_min,
                    dens_mean,
                    dens_max,
                    &r, &g, &b);
      eps_file << i - 0.5*sq_overlap << " " << j - 0.5*sq_overlap  << " sq ";
      eps_file << r << " " << g << " " << b << " srgb f\n";
    }
  }
  write_polygons_to_eps(eps_file, false, false, inset_state);
  eps_file << "showpage\n";
  eps_file << "%%EOF\n";
  eps_file.close();
  return;
}
