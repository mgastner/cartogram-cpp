#include "constants.hpp"
#include "inset_state.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/srs/projections/proj4.hpp>
#include <boost/geometry/srs/transformation.hpp>

namespace bg = boost::geometry;
namespace srs = bg::srs;

class AlbersProjector
{
  using trans_t = srs::transformation<>;

public:
  AlbersProjector(double lambda0, double phi0, double phi1, double phi2)
      : lambda0_(lambda0), phi1_(phi1),
        cylindrical_(std::abs(phi1 + phi2) < 1e-6)
  {
    if (cylindrical_)
      return;

    // Write up parameters to be used for Albers projection
    // NOTE: R (radius) is set to 1 for our conversion
    const std::string src = "+proj=longlat +R=1 +no_defs";
    const std::string dst =
      "+proj=aea +lat_1=" + std::to_string(phi1) +
      " +lat_2=" + std::to_string(phi2) + " +lat_0=" + std::to_string(phi0) +
      " +lon_0=" + std::to_string(lambda0) + " +R=1 +no_defs";

    trans_ = std::make_unique<trans_t>(srs::proj4(src), srs::proj4(dst));
  }

  Point operator()(const Point &p) const
  {
    if (cylindrical_) {
      // If n = 0 (i.e., phi_1 = -phi_2), the Albers projection becomes a
      // cylindrical equal-area projection with standard parallel phi_1. The
      // formula is at:
      // https://en.wikipedia.org/wiki/Cylindrical_equal-area_projection
      const double lon = p.x() * pi / 180.0;
      const double lat = p.y() * pi / 180.0;
      const double x = (lon - lambda0_) * std::cos(phi1_);
      const double y = std::sin(lat) / std::cos(phi1_);
      return {x, y};
    }

    bg::model::point<double, 2, bg::cs::geographic<bg::degree>> in(
      p.x(),
      p.y());
    bg::model::point<double, 2, bg::cs::cartesian> out;

    (*trans_).forward(in, out);

    return {bg::get<0>(out), bg::get<1>(out)};
  }

private:
  double lambda0_;
  double phi1_;
  bool cylindrical_;
  std::unique_ptr<trans_t> trans_;
};

void InsetState::adjust_for_dual_hemisphere()
{
  // Determine the maximum longitude in the western hemisphere and the minimum
  // longitude in the eastern hemisphere
  double max_lon_west = -dbl_inf;
  double min_lon_east = dbl_inf;
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const auto bb = pwh.bbox();
      const double xmax = bb.xmax();
      const double xmin = bb.xmin();
      max_lon_west = xmax < 0 ? std::max(xmax, max_lon_west) : max_lon_west;
      min_lon_east = xmin >= 0 ? std::min(xmin, min_lon_east) : min_lon_east;
    }
  }

  // Set transformation (translation) values to +360 for longitude
  Transformation translate(CGAL::TRANSLATION, CGAL::Vector_2<Scd>(360, 0));

  // - If min_lon_east == max_lon_west, the whole inset is contained in either
  //   only the western or only the eastern hemisphere
  // - If max_lon_west < -180.0, all polygons that are partly in the western
  //   hemisphere also are partly in the eastern hemisphere
  // - If min_lon_east > 180.0, all polygons that are partly in the eastern
  //   hemisphere also are partly in the western hemisphere
  // - If min_lon_east - max_lon_west < 180, the inset cannot fit in 1
  //   hemisphere
  if (
    max_lon_west >= -180.0 && min_lon_east <= 180.0 &&
    min_lon_east - max_lon_west >= 180) {

    // Iterate over GeoDivs
    for (auto &gd : geo_divs_) {

      // Iterate over Polygon_with_holes
      for (auto &pwh : gd.ref_to_polygons_with_holes()) {
        auto &outer_boundary = pwh.outer_boundary();

        // If pwh is in the western hemisphere
        if (pwh.bbox().xmin() < 0) {
          outer_boundary = transform(translate, outer_boundary);

          // Iterate over holes
          for (auto &h : pwh.holes()) {
            h = transform(translate, h);
          }
        }
      }
    }
  }
}

void InsetState::apply_albers_projection()
{
  // Adjust the longitude coordinates if the inset spans both the eastern and
  // western hemispheres
  adjust_for_dual_hemisphere();

  // Recalculate the bbox after dual hemisphere adjustment
  const auto bb = bbox();

  // Declarations for albers_formula()
  // const double min_lon = (bb.xmin() * pi) / 180;
  // const double min_lat = (bb.ymin() * pi) / 180;
  // const double max_lon = (bb.xmax() * pi) / 180;
  // const double max_lat = (bb.ymax() * pi) / 180;
  const double min_lon = bb.xmin();
  const double min_lat = bb.ymin();
  const double max_lon = bb.xmax();
  const double max_lat = bb.ymax();

  // Reference longitude and latitude
  const double lambda_0 = 0.5 * (min_lon + max_lon);
  const double phi_0 = 0.5 * (min_lat + max_lat);

  // Standard parallels
  const double phi_1 = 0.5 * (phi_0 + max_lat);
  const double phi_2 = 0.5 * (phi_0 + min_lat);

  static thread_local AlbersProjector proj(lambda_0, phi_0, phi_1, phi_2);

  transform_points([&](const Point &q) {
    return proj(q);
  });
}
