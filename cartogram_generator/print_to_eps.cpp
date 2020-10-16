#include "map_state.h"
#include <fstream>

void print_to_eps(std::string eps_name, MapState *map_state)
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

  eps_file << "100 100 moveto 350 0 rlineto 0 350 rlineto -350 0 rlineto "
           << "closepath stroke\n";
  eps_file << "showpage\n";
  eps_file << "%%EOF\n";
  eps_file.close();
  return;
}
