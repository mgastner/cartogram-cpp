#include "inset_state.h"
#include "constants.h"
#include <cmath>

// Functions to project map with the Smyth equal-surface projection (also   
// known as Craster rectangular projection). See                            
// https://en.wikipedia.org/wiki/Cylindrical_equal-area_projection          
// The purpose is to create a projection that has a 2:1 aspect ratio so that
// the Fourier transforms work optimally when padding is reduced to zero.   

double project_x_to_smyth(double x)
{
  return x * std::sqrt(2.0 * pi) / 180.0;
}

double project_y_to_smyth(double y)
{
  return std::sin(y * pi / 180.0) * std::sqrt(0.5 * pi);
}

void project_to_smyth_equal_surface(InsetState *inset_state)
{

  std::vector<GeoDiv> new_geo_divs;

  for (auto gd : inset_state->geo_divs()) {

    // For each GeoDiv
    GeoDiv new_gd(gd.id());

    for (auto pwh : gd.polygons_with_holes()) {
      // For each polygon with holes

      Polygon old_ext_ring = pwh.outer_boundary();
      Polygon new_ext_ring;

      for (unsigned int i = 0; i < old_ext_ring.size(); i++) {

        // Update exterior ring coordinates
        double new_ext_ring_point_x =
          project_x_to_smyth(old_ext_ring[i][0]);

        double new_ext_ring_point_y =
          project_y_to_smyth(old_ext_ring[i][1]);

        new_ext_ring.push_back(Point(new_ext_ring_point_x,
                                     new_ext_ring_point_y));
      }
      std::vector<Polygon> hole_v;
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); hci++) {
        Polygon old_hole = *hci;
        Polygon new_hole;
        for (unsigned int i = 0; i < old_hole.size(); i++) {

          // Update hole coordinates
          double new_hole_point_x =
            project_x_to_smyth(old_hole[i][0]);

          double new_hole_point_y =
            project_y_to_smyth(old_hole[i][1]);
            
          new_hole.push_back(Point(new_hole_point_x,
                                   new_hole_point_y));
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
  inset_state->set_geo_divs(new_geo_divs);

  return;
}

/* Function for the projection from Smyth-Craster coordinates to longitude/  */
/* latitude. We assume that the Smyth-Craster coordinates have been scaled   */
/* to fit in the box [0, lx] * [0, ly].                                      */

double project_x_from_smyth(double x, int lx)
{
  return x * 360.0 / lx - 180.0;
}

double project_y_from_smyth(double y, int ly)
{
  return 180.0 * std::asin((2.0 * y / ly) - 1) / pi;
}

void project_from_smyth_equal_surface(InsetState *inset_state)
{

  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();

  std::vector<GeoDiv> new_geo_divs;

  for (auto gd : inset_state->geo_divs()) {

    // For each GeoDiv
    GeoDiv new_gd(gd.id());

    for (auto pwh : gd.polygons_with_holes()) {
      // For each polygon with holes

      Polygon old_ext_ring = pwh.outer_boundary();
      Polygon new_ext_ring;

      for (unsigned int i = 0; i < old_ext_ring.size(); i++) {

        // Update exterior ring coordinates
        double new_ext_ring_point_x =
          project_x_from_smyth(old_ext_ring[i][0], lx);

        double new_ext_ring_point_y =
          project_y_from_smyth(old_ext_ring[i][1], ly);

        new_ext_ring.push_back(Point(new_ext_ring_point_x,
                                     new_ext_ring_point_y));
      }
      std::vector<Polygon> hole_v;
      for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); hci++) {
        Polygon old_hole = *hci;
        Polygon new_hole;
        for (unsigned int i = 0; i < old_hole.size(); i++) {

          // Update hole coordinates
          double new_hole_point_x =
            project_x_from_smyth(old_hole[i][0], lx);

          double new_hole_point_y =
            project_y_from_smyth(old_hole[i][1], ly);
            
          new_hole.push_back(Point(new_hole_point_x,
                                   new_hole_point_y));
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
  inset_state->set_geo_divs(new_geo_divs);

  return;
}
