#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <iostream>
#include <fstream>
#include <vector>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_with_holes_2<K>                       PolygonWH;
typedef CGAL::Polygon_2<K>                                  Polygon_2;
typedef CGAL::Point_2<K>                                    Point;

class GeoDiv {

  private:
    std::string id;
    std::vector<PolygonWH> polygons_wh;

  public:
    GeoDiv(std::string i);
    std::vector<PolygonWH> GetPolygons();
    void AddPolygon(PolygonWH pgn_wh);

};

GeoDiv::GeoDiv(std::std::string i) { // constructor
    id = i;
}

std::vector<PolygonWH> GeoDiv::GetPolygons() { // return polygons
  return polygons_wh;
}

void GeoDiv::AddPolygon(PolygonWH pgn_wh) { // add a polygon
  polygons_wh.push_back(pgn_wh);
}

std::string HexToRGB (std::string hex)
{

  if (hex[0] == '#') {
    hex.erase(0, 1);
  }

  int hex_int = stoi(hex, nullptr, 16);

  double red = ((hex_int >> 16) & 0xFF) / 255.0; // Extract the RR byte
  double green = ((hex_int >> 8) & 0xFF) / 255.0; // Extract the GG byte
  double blue = ((hex_int) & 0xFF) / 255.0; // Extract the BB byte

  // from https://gist.github.com/bert/998020

  std::string rgb;

  rgb.append(to_string(red));
  rgb.append(" ");
  rgb.append(to_string(green));
  rgb.append(" ");
  rgb.append(to_string(blue));
  rgb.append(" ");

  return rgb;

}

void ps_figure (std::string ps_name, std::vector<GeoDiv> map)
{
  cout << "Entered ps_figure()" << endl;
  ofstream ps_file(ps_name);

  /**************************** Postscript header. ***************************/

  ps_file << "%!PS-Adobe-2.0 EPSF-2.0" << endl;
  ps_file << "%%Title: " << ps_name << endl;
  ps_file << "%%Creator: Michael T. Gastner et al." << endl;
  ps_file << "%%BoundingBox: 0 0 20 20" << endl;
  ps_file << "%%Magnification: 1.0000" << endl;
  ps_file << "%%EndComments" << endl;
  ps_file << "/m {moveto} def" << endl;
  ps_file << "/l {lineto} def" << endl;
  ps_file << "/s {stroke} def" << endl;
  ps_file << "/n {newpath} def" << endl;
  ps_file << "/c {closepath} def" << endl;
  ps_file << "/f {fill} def" << endl;
  ps_file << "/SLW {setlinewidth} def" << endl;
  ps_file << "/SGRY {setgray} def" << endl;
  ps_file << "/SRGB {setrgbcolor} def" << endl;

  /**************************** Plot the regions. ****************************/

  // for each geodiv within map
  for (GeoDiv region : map) {

    // for each polygon within geodiv
    for (PolygonWH pgnWH : region.GetPolygons()) {

    // setting line width and gray scale to black
    ps_file << "0 SGRY" << endl;
    ps_file << "0.1 SLW" << endl;


      /****** PLOTTING POLYGON (OUTER BOUNDARY) ******/

      Polygon_2 pgn = pgnWH.outer_boundary();

      if (pgn.size() > 1) {

        // starting a new line
        ps_file << "n" << endl;

        if (pgn.is_clockwise_oriented()) {
          //POLYGON CLOCKWISE - MUST SWITCH ORIENTATION***//

          // moving to starting co-ordinates
          ps_file << pgn[pgn.size() - 1][0] << " " << pgn[pgn.size() - 1][1] << " m" << endl;

          // plotting each point in polygon
          for (int i = pgn.size() - 1; i >= 0; --i) {
            ps_file << pgn[i][0] << " " << pgn[i][1] << " l" << endl;

          }
        } else {
          //***HOLES COUNTER-CLOCKWISE - CAN PRINT AS IS***//

          // moving to starting co-ordinates
          ps_file << pgn[0][0] << " " << pgn[0][1] << " m" << endl;

          // plotting each point in polygon
          for (int i = 1; i < pgn.size(); ++i) {
            ps_file << pgn[i][0] << " " << pgn[i][1] << " l" << endl;

          }
        }

        // closing path
        ps_file << "c" << endl;

        /*********** PLOTTING HOLES ***********/

        std::vector<Polygon_2> holes(pgnWH.holes_begin(), pgnWH.holes_end());

        for (Polygon_2 hole : holes) {

          if (hole.size() > 0) {

            // DO NOT START NEW LINE HERE

            if (hole.is_clockwise_oriented()) {

              //***HOLES CLOCKWISE - CAN PRINT AS IS***//

              ps_file << hole[0][0] << " " << hole[0][1] << " m" << endl;

              // plotting points in hole
              for (int i = 1; i < hole.size(); ++i) {
                ps_file << hole[i][0] << " " << hole[i][1] << " l" << endl;
              }

            } else {
              //***HOLE COUNTER-CLOCKWISE - MUST SWITCH ORIENTATION***//

              // printing from opposite direction
              ps_file << hole[hole.size() - 1][0] << " " << hole[hole.size() - 1][1] << " m" << endl;

              // plotting points in hole in opposite direction
              for (int i = hole.size() - 1; i >= 0; --i) {
                ps_file << hole[i][0] << " " << hole[i][1] << " l" << endl;
              }
            }

            // closing path
            ps_file << "c" << endl;

          }

        }

        // saving path before filling it
        ps_file << "gsave" << endl;

        // filling path
        ps_file << HexToRGB("#34495e") << "SRGB f" << endl;

        // restoring path and stroking outline
        ps_file << "grestore" << endl;
        ps_file << "s" << endl;

      }

    }

  }

  ps_file << "showpage" << endl;
  ps_file.close();

  return;
}

int main()
{
  // create a polygon with three holes
  Polygon_2 outer_polygon;
  // clockwise polygon
  outer_polygon.push_back(Point(0,0)); outer_polygon.push_back(Point(0,9));
  outer_polygon.push_back(Point(9, 9)); outer_polygon.push_back(Point(9,0));
  std::vector<Polygon_2> holes(2);

  // counter-clockwise points
  holes[0].push_back(Point(6,6)); holes[0].push_back(Point(6,8));
  holes[0].push_back(Point(8, 8)); holes[0].push_back(Point(8,6));

  // clockwise points
  holes[1].push_back(Point(2,4)); holes[1].push_back(Point(2,2));
  holes[1].push_back(Point(4,2)); holes[1].push_back(Point(4,4));

  PolygonWH pgnWH(outer_polygon, holes.begin(), holes.end());

  GeoDiv region("India");
  region.AddPolygon(pgnWH);

  std::vector<GeoDiv> map;
  map.push_back(region);

  CGAL::set_pretty_mode(cout);
  cout << map[0].GetPolygons()[0];

  // And draw it.
  ps_figure("test.eps", map);

  return EXIT_SUCCESS;
}
