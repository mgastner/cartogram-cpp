#include "inset_state.h"
#include "constants.h"
#include <cmath>

// Functions to project map with the Smyth equal-surface projection (also
// known as Craster rectangular projection):
// https://en.wikipedia.org/wiki/Cylindrical_equal-area_projection
// The purpose is to create a projection that has a 2:1 aspect ratio so that
// the Fourier transforms work optimally when padding is reduced to zero.
double project_x_to_smyth(const double x)
{
  return x * std::sqrt(2.0 * pi) / 180.0;
}

double project_y_to_smyth(const double y)
{
  return std::sin(y * pi / 180.0) * std::sqrt(0.5 * pi);
}

void project_to_smyth_equal_surface(InsetState *inset_state)
{
  std::vector<GeoDiv> new_geo_divs;
  for (const auto &gd : inset_state->geo_divs()) {
    GeoDiv new_gd(gd.id());
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto old_ext_ring = pwh.outer_boundary();
      Polygon new_ext_ring;
      for (unsigned int i = 0; i < old_ext_ring.size(); ++i) {

        // Update exterior ring coordinates
        const double new_ext_ring_point_x =
          project_x_to_smyth(old_ext_ring[i].x());
        const double new_ext_ring_point_y =
          project_y_to_smyth(old_ext_ring[i].y());
        new_ext_ring.push_back(Point(new_ext_ring_point_x,
                                     new_ext_ring_point_y));
      }
      std::vector<Polygon> hole_v;
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        Polygon new_hole;
        for (unsigned int i = 0; i < h->size(); ++i) {

          // Update hole coordinates
          const double new_hole_point_x =
            project_x_to_smyth((*h)[i].x());
          const double new_hole_point_y =
            project_y_to_smyth((*h)[i].y());
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

/* Functions for the projection from Smyth-Craster coordinates to longitude/ */
/* latitude. We assume that the Smyth-Craster coordinates have been scaled   */
/* to fit in the box [0, lx] * [0, ly].                                      */
double project_x_from_smyth(const double x, const int lx)
{
  return (x * 360.0 / lx) - 180.0;
}

double project_y_from_smyth(const double y, const int ly)
{
  return 180.0 * std::asin((2.0 * y / ly) - 1) / pi;
}

void project_from_smyth_equal_surface(InsetState *inset_state)
{
  const unsigned int lx = inset_state->lx();
  const unsigned int ly = inset_state->ly();
  std::vector<GeoDiv> new_geo_divs;
  for (const auto &gd : inset_state->geo_divs()) {
    GeoDiv new_gd(gd.id());
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto old_ext_ring = pwh.outer_boundary();
      Polygon new_ext_ring;
      for (unsigned int i = 0; i < old_ext_ring.size(); ++i) {

        // Update exterior ring coordinates
        double new_ext_ring_point_x =
          project_x_from_smyth(old_ext_ring[i].x(), lx);
        double new_ext_ring_point_y =
          project_y_from_smyth(old_ext_ring[i].y(), ly);
        new_ext_ring.push_back(Point(new_ext_ring_point_x,
                                     new_ext_ring_point_y));
      }
      std::vector<Polygon> hole_v;
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        Polygon new_hole;
        for (unsigned int i = 0; i < h->size(); ++i) {

          // Update hole coordinates
          const double new_hole_point_x =
            project_x_from_smyth((*h)[i].x(), lx);
          const double new_hole_point_y =
            project_y_from_smyth((*h)[i].y(), ly);
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
