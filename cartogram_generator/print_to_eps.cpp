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
  eps_file << "%!PS-Adobe-3.0 EPSF-3.0 " << std::endl;
  eps_file << "%!PS-Adobe-2.0 EPSF-2.0" << std::endl;
  eps_file << "%%Title: " << eps_name << std::endl;
  eps_file << "%%Creator: Michael T. Gastner et al." << std::endl;
  eps_file << "%%For: Humanity" << std::endl;
  eps_file << "%%Copyright: License CC BY-NC-ND 2.0" << std::endl;
  eps_file << "%%CreationDate: " << asctime(ti);
  eps_file << "%%Pages: 1" << std::endl;
  eps_file << "%%BoundingBox: 0 0 "
           << map_state->get_lx() << " "
           << map_state->get_ly() << std::endl;
  eps_file << "%%Magnification: 1.0000" << std::endl;
  eps_file << "%%EndComments" << std::endl;

  eps_file << "100 100 moveto 350 0 rlineto 0 350 rlineto -350 0 rlineto "
           << "closepath stroke" << std::endl;

  eps_file << "showpage" << std::endl;
  eps_file << "%%EOF" << std::endl;
  eps_file.close();
  return;
}
