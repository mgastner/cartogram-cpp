#include "densify.h"
#include "densification_points.h"

std::vector<GeoDiv> densified_geo_divs(std::vector<GeoDiv> gds)
{
  std::cerr << "Densifying" << std::endl;
  std::vector<GeoDiv> dens_gds;
  for (const auto &gd : gds) {
    GeoDiv dens_gd(gd.id());
    for (const auto &pwh : gd.polygons_with_holes()) {
      const Polygon outer = pwh.outer_boundary();
      Polygon dens_outer;

      // Iterate over each point in the outer boundary of the polygon
      for (size_t i = 0; i < outer.size(); ++i) {

        // The segment defined by points `a` and `b` is to be densified.
        // `b` should be the point immediately after `a`, unless `a` is the
        // final point of the boundary, in which case `b` should be the first
        // point.
        Point a = outer[i];
        Point b = (i == outer.size() - 1) ? outer[0] : outer[i + 1];

        // Densify the segment
        std::vector<Point> dens_outer_pts = densification_points(a, b);

        // Push all points. Omit the last point because it will be included
        // in the next iteration. Otherwise, we would have duplicated points
        // in the polygon.
        for (size_t j = 0; j < (dens_outer_pts.size() - 1); ++j) {
          dens_outer.push_back(dens_outer_pts[j]);
        }
      }
      std::vector<Polygon> dens_holes;
      std::vector<Polygon> holes(pwh.holes_begin(), pwh.holes_end());
      for (Polygon hole : holes) {
        Polygon dens_hole;

        // Iterate over each point in the hole
        for (size_t i = 0; i < hole.size(); ++i) {

          // `c` and `d` are determined in the same way as `a` and `b` above
          Point c = hole[i];
          Point d = (i == hole.size() - 1) ? hole[0] : hole[i + 1];
          std::vector<Point> dens_hole_pts = densification_points(c, d);
          for (size_t j = 0; j < (dens_hole_pts.size() - 1); ++j) {
            dens_hole.push_back(dens_hole_pts[j]);
          }
        }
        dens_holes.push_back(dens_hole);
      }
      Polygon_with_holes dens_pwh(dens_outer,
                                  dens_holes.begin(),
                                  dens_holes.end());
      dens_gd.push_back(dens_pwh);
    }
    dens_gds.push_back(dens_gd);
  }
  return dens_gds;
}
