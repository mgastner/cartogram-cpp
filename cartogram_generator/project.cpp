#include <boost/multi_array.hpp>
#include <iostream>
#include <vector>
#include "interpolate_bilinearly.h"

#include "project.h"

void project(MapState *map_state)
{
  const unsigned int lx = map_state->lx();
  const unsigned int ly = map_state->ly();
  boost::multi_array<XYPoint, 2> &proj = *map_state->proj();

  // Calculate displacement from proj array
  boost::multi_array<double, 2> xdisp(boost::extents[lx][ly]);
  boost::multi_array<double, 2> ydisp(boost::extents[lx][ly]);
  for (unsigned int i = 0; i < lx; i++) {
    for (unsigned int j=0; j<ly; j++) {
      xdisp[i][j] = proj[i][j].x - i - 0.5;
      ydisp[i][j] = proj[i][j].y - j - 0.5;
    }
  }
  std::vector<GeoDiv> new_geo_divs;
  for (auto gd : map_state->geo_divs()) {

    // For each GeoDiv
    GeoDiv new_gd(gd.id());

    for (auto pwh : gd.polygons_with_holes()) {
      // For each polygon with holes

      Polygon old_ext_ring = pwh.outer_boundary();
      Polygon new_ext_ring;

      for (unsigned int i = 0; i < old_ext_ring.size(); i++) {

        // Update exterior ring coordinates
        double old_ext_intp_x =
          interpolate_bilinearly(old_ext_ring[i][0], old_ext_ring[i][1],
                                 &xdisp, 'x',
                                 lx, ly);
        double old_ext_intp_y =
          interpolate_bilinearly(old_ext_ring[i][0], old_ext_ring[i][1],
                                 &ydisp, 'y',
                                 lx, ly);
        new_ext_ring.push_back(Point(old_ext_intp_x +old_ext_ring[i][0],
                                     old_ext_intp_y + old_ext_ring[i][1]));
      }
      std::vector<Polygon> hole_v;
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); hci++) {
        Polygon old_hole = *hci;
        Polygon new_hole;
        for (unsigned int i = 0; i < old_hole.size(); i++) {

          // Update hole coordinates
          double old_hole_intp_x =
            interpolate_bilinearly(old_hole[i][0], old_hole[i][1],
                                   &xdisp, 'x',
                                   lx, ly);
          double old_hole_intp_y =
            interpolate_bilinearly(old_hole[i][0], old_hole[i][1],
                                   &ydisp, 'y', lx, ly);
            new_hole.push_back(Point(old_hole_intp_x + old_hole[i][0],
                                     old_hole_intp_y + old_hole[i][1]));
        }
        hole_v.push_back(new_hole);
      }
      const Polygon_with_holes new_pwh(new_ext_ring,
                                       hole_v.begin(),
                                       hole_v.end());
      new_gd.push_back(new_pwh);
    }
    new_geo_divs.push_back(new_gd);
  }
  map_state->set_geo_divs(new_geo_divs);
  return;
}

std::vector<int> find_graticule_point(const int x,
                                      const int y,
                                      const int lx,
                                      const int ly)
{
  if (x < 0 || x > lx || y < 0 || y > ly) {
    std::cout << "ERROR: coordinate outside bounding box in find_graticule_point()." << "\n";
    std::cout << "x=" << x << ", y=" << y << "\n";
    exit(1);
  }

  std::vector<int> graticule_point;

  graticule_point.push_back(int(floor(x)));
  graticule_point.push_back(int(floor(y)));

  return graticule_point;
}

void project_graticule(MapState *map_state)
{
  const unsigned int lx = map_state->lx();
  const unsigned int ly = map_state->ly();
  boost::multi_array<XYPoint, 2> &proj = *map_state->proj();
  
  boost::multi_array<XYPoint, 2> &graticule_points = *map_state->graticule_points();
  
  // Resize multi array if running for the first time
  if (map_state->n_finished_integrations() == 0) {
    graticule_points.resize(boost::extents[lx + 1][ly + 1]);
  }

  for (unsigned int i = 0; i < lx + 1; i++){
    for (unsigned int j = 0; j < ly + 1; j++){
      graticule_points[i][j].x = i;
      graticule_points[i][j].y = j;
    }
  }

  // Calculate displacement from proj array
  boost::multi_array<double, 2> xdisp(boost::extents[lx][ly]);
  boost::multi_array<double, 2> ydisp(boost::extents[lx][ly]);
  for (unsigned int i = 0; i < lx; i++) {
    for (unsigned int j=0; j<ly; j++) {
      xdisp[i][j] = proj[i][j].x - i - 0.5;
      ydisp[i][j] = proj[i][j].y - j - 0.5;
    }
  }

  // Project graticules
  for (unsigned int i = 0; i < lx; i++){
    for (unsigned int j = 0; j < ly; j++){
      double new_x =
        interpolate_bilinearly(graticule_points[i][j].x,
                               graticule_points[i][j].y,
                               &xdisp, 'x', lx, ly);
      double new_y =
        interpolate_bilinearly(graticule_points[i][j].x,
                               graticule_points[i][j].y,
                               &ydisp, 'y', lx, ly);
      graticule_points[i][j].x = new_x;
      graticule_points[i][j].y = new_y;
    }
  }
  return;
}