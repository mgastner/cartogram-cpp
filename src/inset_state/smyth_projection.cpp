#include "../constants.h"
#include "../inset_state.h"
#include <math.h>

// Functions to project map with the Smyth equal-surface projection (also
// known as Craster rectangular projection):
// https://en.wikipedia.org/wiki/Cylindrical_equal-area_projection
// The purpose is to create a projection that has a 2:1 aspect ratio so that
// the Fourier transforms work optimally when padding is reduced to zero.
Point apply_smyth_craster_projection_on_point(Point p1)
{
  return Point(p1.x() * std::sqrt(2.0 * pi) / 180.0,
               std::sin(p1.y() * pi / 180.0) * std::sqrt(0.5 * pi));
}

void InsetState::apply_smyth_craster_projection()
{
  // Iterate over GeoDivs
  for (auto &gd : geo_divs_) {

    // Iterate over Polygon_with_holes
    for (auto &pwh : *gd.ref_to_polygons_with_holes()) {

      // Get outer boundary
      auto &outer_boundary = *(&pwh.outer_boundary());

      // Iterate over outer boundary's coordinates
      for (auto &coords_outer : outer_boundary) {

        // Assign outer boundary's coordinates to transformed coordinates
        coords_outer = apply_smyth_craster_projection_on_point(coords_outer);
      }

      // Iterate over holes
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {

        // Iterate over hole's coordinates
        for (auto &coords_hole : *h) {

          // Assign hole's coordinates to transformed coordinates
          coords_hole = apply_smyth_craster_projection_on_point(coords_hole);
        }
      }
    }
  }
  return;
}

// Functions for the projection from Smyth-Craster coordinates to longitude
// latitude. We assume that the Smyth-Craster coordinates have been scaled
// to fit in the box [0, lx] * [0, ly].
Point reverse_smyth_craster_projection_on_point(Point p1,
                                                 const unsigned int lx,
                                                 const unsigned int ly)
{
  return Point((p1.x() * 360.0 / lx) - 180.0,
               180.0 * std::asin((2.0 * p1.y() / ly) - 1) / pi);
}

void InsetState::revert_smyth_craster_projection()
{
  // Iterate over GeoDivs
  for (auto &gd : geo_divs_) {

    // Iterate over Polygon_with_holes
    for (auto &pwh : *gd.ref_to_polygons_with_holes()) {

      // Get outer boundary
      auto &outer_boundary = *(&pwh.outer_boundary());

      // Iterate over outer boundary's coordinates
      for (auto &coords_outer : outer_boundary) {

        // Assign outer boundary's coordinates to transformed coordinates
        coords_outer = reverse_smyth_craster_projection_on_point(coords_outer,
                                                                 lx_,
                                                                 ly_);
      }

      // Iterate over holes
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {

        // Iterate over hole's coordinates
        for (auto &coords_hole : *h) {

          // Assign hole's coordinates to transformed coordinates
          coords_hole = reverse_smyth_craster_projection_on_point(coords_hole,
                                                                  lx_,
                                                                  ly_);
        }
      }
    }
  }
  return;
}
