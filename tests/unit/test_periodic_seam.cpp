#define BOOST_TEST_MODULE test_periodic_seam
#include "constants.hpp"
#include "inset_state.hpp"
#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <set>
#include <vector>

// With --periodic_x, regions can flow across the seam at x = 0 and x = lx.
// The tests check the treatment of the seam outside the flow: the output
// cuts polygons at the sides of the frame, simplification keeps the vertices
// where the input cuts regions at the antimeridian, and an input ring that
// crosses the antimeridian without being cut is found.

namespace
{

Polygon ring_from(const std::vector<std::pair<double, double>> &xy)
{
  Polygon p;
  for (const auto &[x, y] : xy) {
    p.push_back(Point(x, y));
  }
  return p;
}

double area_of(const GeoDiv &gd)
{
  double a = 0.0;
  for (const auto &pwh : gd.polygons_with_holes()) {
    a += pwh.outer_boundary().area();
    for (const auto &h : pwh.holes()) {
      a += h.area();
    }
  }
  return a;
}

// The y-coordinates of the vertices of gd on the line x = x0
std::multiset<double> ys_on_line(const GeoDiv &gd, const double x0)
{
  std::multiset<double> ys;
  auto add = [&](const Polygon &ring) {
    for (const Point &p : ring) {
      if (!(p.x() < x0) && !(p.x() > x0)) {
        ys.insert(p.y());
      }
    }
  };
  for (const auto &pwh : gd.polygons_with_holes()) {
    add(pwh.outer_boundary());
    for (const auto &h : pwh.holes()) {
      add(h);
    }
  }
  return ys;
}

}  // namespace

// A polygon with a hole crosses x = lx, a second polygon crosses x = 0, and a
// third lies inside the frame. After wrap_into_frame(), every vertex lies in
// the frame, each GeoDiv keeps its area, the inner polygon is unchanged, and
// the pieces on the two sides of a cut have the same cut points.
BOOST_AUTO_TEST_CASE(polygons_across_the_seam_are_cut_into_the_frame)
{
  constexpr unsigned int lx = 32;
  Arguments args{};
  args.world = true;
  args.periodic_x = true;
  InsetState inset("W", args);

  const Polygon hole =
    ring_from({{30.0, 6.0}, {30.0, 10.0}, {36.5, 10.0}, {36.5, 6.0}});
  const std::vector<Polygon> holes = {hole};
  GeoDiv a("A");
  a.push_back(Polygon_with_holes(
    ring_from({{28.0, 4.0}, {38.25, 4.0}, {38.25, 12.0}, {28.0, 12.0}}),
    holes.begin(),
    holes.end()));
  GeoDiv b("B");
  b.push_back(Polygon_with_holes(
    ring_from({{-5.5, 2.0}, {3.0, 2.25}, {3.0, 5.0}, {-5.5, 5.0}})));
  GeoDiv c("C");
  const Polygon inner =
    ring_from({{10.0, 10.0}, {12.0, 10.0}, {12.0, 12.0}, {10.0, 12.0}});
  c.push_back(Polygon_with_holes(inner));
  for (const GeoDiv &gd : {a, b, c}) {
    inset.push_back(gd);
  }
  inset.set_grid_dimensions(lx, 16);
  const std::vector<double> areas_before = {
    area_of(a),
    area_of(b),
    area_of(c)};

  inset.wrap_into_frame();

  const auto &gds = inset.geo_divs();
  BOOST_TEST_REQUIRE(gds.size() == 3u);
  for (std::size_t i = 0; i < gds.size(); ++i) {
    BOOST_TEST(
      std::fabs(area_of(gds[i]) - areas_before[i]) < 1e-12 * areas_before[i]);
    for (const auto &pwh : gds[i].polygons_with_holes()) {
      for (const Point &p : pwh.outer_boundary()) {
        BOOST_TEST(!(p.x() < 0.0));
        BOOST_TEST(!(p.x() > lx));
      }
    }
  }

  // A and B are cut into two pieces each; C is unchanged
  BOOST_TEST(gds[0].polygons_with_holes().size() == 2u);
  BOOST_TEST(gds[1].polygons_with_holes().size() == 2u);
  BOOST_TEST_REQUIRE(gds[2].polygons_with_holes().size() == 1u);
  const Polygon &inner_after =
    gds[2].polygons_with_holes()[0].outer_boundary();
  BOOST_TEST(
    std::equal(
      inner_after.begin(),
      inner_after.end(),
      inner.begin(),
      inner.end()));

  // The cut points on x = lx and on x = 0 are the same points of the seam
  for (std::size_t i = 0; i < 2; ++i) {
    const auto west = ys_on_line(gds[i], 0.0);
    const auto east = ys_on_line(gds[i], lx);
    BOOST_TEST(!west.empty());
    BOOST_TEST((west == east));
  }
}

// Two pieces of a region, cut by the input at the antimeridian, with 41
// vertices each on the cut line. After the projection and the rescaling, the
// cut vertices lie exactly on x = 0 and x = lx, and a strong simplification
// keeps all of them, so that the two copies of the cut line stay identical.
BOOST_AUTO_TEST_CASE(vertices_on_the_antimeridian_survive_simplification)
{
  Arguments args{};
  args.world = true;
  args.periodic_x = true;
  args.n_grid_rows_or_cols = 64;
  InsetState inset("W", args);

  std::vector<std::pair<double, double>> east, west;
  for (int k = 0; k <= 40; ++k) {
    const double lat = -30.0 + 1.5 * k;
    east.emplace_back(180.0, lat);
    west.emplace_back(-180.0, lat);
  }
  for (int k = 40; k >= 0; --k) {
    const double lat = -30.0 + 1.5 * k;
    east.emplace_back(170.0 - 0.01 * (k % 2), lat);
    west.emplace_back(-170.0 + 0.01 * (k % 2), lat);
  }
  Polygon e = ring_from(east);
  Polygon w = ring_from(west);
  for (Polygon *p : {&e, &w}) {
    if (p->is_clockwise_oriented()) {
      p->reverse_orientation();
    }
  }
  GeoDiv gd("R");
  gd.push_back(Polygon_with_holes(e));
  gd.push_back(Polygon_with_holes(w));
  inset.push_back(gd);

  inset.apply_smyth_craster_projection();
  inset.rescale_map();
  const double dlx = inset.lx();
  const auto west_before = ys_on_line(inset.geo_divs()[0], 0.0);
  const auto east_before = ys_on_line(inset.geo_divs()[0], dlx);
  BOOST_TEST_REQUIRE(west_before.size() == 41u);
  BOOST_TEST_REQUIRE((west_before == east_before));

  inset.simplify(10);

  BOOST_TEST(inset.n_points() < 164u);
  BOOST_TEST((ys_on_line(inset.geo_divs()[0], 0.0) == west_before));
  BOOST_TEST((ys_on_line(inset.geo_divs()[0], dlx) == east_before));
}

// An edge from 170 to -170 degrees crosses the antimeridian; the same region
// cut at 180 degrees does not
BOOST_AUTO_TEST_CASE(edge_across_the_antimeridian_is_found)
{
  Arguments args{};
  args.world = true;
  args.periodic_x = true;
  InsetState uncut("W", args);
  GeoDiv gd("Fiji");
  gd.push_back(Polygon_with_holes(ring_from(
    {{170.0, -20.0}, {-170.0, -20.0}, {-170.0, -10.0}, {170.0, -10.0}})));
  uncut.push_back(gd);
  BOOST_TEST(uncut.geo_div_with_edge_across_antimeridian() == "Fiji");

  InsetState cut("W", args);
  GeoDiv gd_cut("Fiji");
  gd_cut.push_back(Polygon_with_holes(ring_from(
    {{170.0, -20.0}, {180.0, -20.0}, {180.0, -10.0}, {170.0, -10.0}})));
  gd_cut.push_back(Polygon_with_holes(ring_from(
    {{-180.0, -20.0}, {-170.0, -20.0}, {-170.0, -10.0}, {-180.0, -10.0}})));
  cut.push_back(gd_cut);
  BOOST_TEST(cut.geo_div_with_edge_across_antimeridian().empty());
}
