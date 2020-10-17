#include "map_state.h"
#include <fstream>

void write_eps(std::string eps_name, MapState *map_state)
{
  std::ofstream eps_file(eps_name);

  // Current time
  time_t tt;
  time(&tt);
  struct tm *ti = localtime(&tt);

  // EPS header
  eps_file << "%!PS-Adobe-3.0 EPSF-3.0\n";
  eps_file << "%%Title: " << eps_name << "\n";
  eps_file << "%%Creator: Michael T. Gastner et al.\n";
  eps_file << "%%For: Humanity\n";
  eps_file << "%%Copyright: License CC BY-NC-ND 2.0\n";
  eps_file << "%%CreationDate: " << asctime(ti);
  eps_file << "%%Pages: 1\n";
  eps_file << "%%BoundingBox: 0 0 "
           << map_state->get_lx() << " "
           << map_state->get_ly() << "\n";
  eps_file << "%%Magnification: 1.0000\n";
  eps_file << "%%EndComments\n";

  // Definitions
  eps_file << "/m {moveto} def\n";
  eps_file << "/l {lineto} def\n";
  eps_file << "/s {stroke} def\n";
  eps_file << "/n {newpath} def\n";
  eps_file << "/c {closepath} def\n";
  eps_file << "/f {fill} def\n";
  eps_file << "/slw {setlinewidth} def\n";
  eps_file << "/sgry {setgray} def\n";
  eps_file << "/srgb {setrgbcolor} def\n";

  // Print polygons
  for (auto gd : map_state->get_geo_divs()) {
    for (auto pwh : gd.get_polygons_with_holes()) {
      CGAL::Polygon_2<K> ext_ring = pwh.outer_boundary();

      // Move to starting coordinates
      eps_file << "n " << ext_ring[0][0] << " " << ext_ring[0][1] << " m\n";

      // Plot each point in exterior ring
      for (int i = 1; i < ext_ring.size(); ++i) {
        eps_file << ext_ring[i][0] << " " << ext_ring[i][1] << " l\n";
      }

      // Close path
      eps_file << "c\n";

      // Plot holes
      typedef typename CGAL::Polygon_with_holes_2<K>::Hole_const_iterator HCI;
      for (HCI hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
        CGAL::Polygon_2<K> hole = *hci;
        eps_file << hole[0][0] << " " << hole[0][1] << " m\n";
        for (int i = 1; i < hole.size(); ++i) {
          eps_file << hole[i][0] << " " << hole[i][1] << " l\n";
        }
        eps_file << "c\n";
      }

      // Save path before filling it
      eps_file << "gsave\n";

      // Fill path
      eps_file << "0.96 0.92 0.70 srgb f\n";

      // Restore path and stroke outline
      eps_file << "grestore\n";
      eps_file << "0 sgry s\n";
    }
  }
  eps_file << "showpage\n";
  eps_file << "%%EOF\n";
  eps_file.close();
  return;
}
