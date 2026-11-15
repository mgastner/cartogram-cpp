#define BOOST_TEST_MODULE test_fill_with_density_lattice_lines
#include "inset_state.hpp"
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// The function fill_with_density() gives each lattice cell the average of the
// target densities of the geographic divisions and of the ocean that overlap
// the cell, weighted by overlap area times area error. The polygons below
// have edges on lattice lines and vertices on lattice points, so that the
// boundaries of the polygons run along the sides of cells or through their
// corners. We compare the density of every cell with the density computed
// from the exact overlap areas. Cells whose only contact with a ring was along
// their sides or at their corners used to get the overlap area 0 with that
// ring, even if they were inside the ring. For example, the cells along the
// left and bottom sides of a rectangle got the ocean density.

namespace
{

constexpr unsigned int lx = 32;
constexpr unsigned int ly = 32;

using Coordinates = std::vector<std::pair<double, double>>;

Polygon ring(const Coordinates &coordinates)
{
  Polygon p;
  for (const auto &[x, y] : coordinates) {
    p.push_back(Point(x, y));
  }
  return p;
}

// Rectangle [x0, x1] x [y0, y1] in counterclockwise order, with a vertex at
// every lattice point on its sides, starting in the middle of the bottom side
Polygon densified_rectangle(
  const int x0,
  const int y0,
  const int x1,
  const int y1)
{
  Polygon p;
  const int x_mid = (x0 + x1) / 2;
  for (int x = x_mid; x < x1; ++x) {
    p.push_back(Point(x, y0));
  }
  for (int y = y0; y < y1; ++y) {
    p.push_back(Point(x1, y));
  }
  for (int x = x1; x > x0; --x) {
    p.push_back(Point(x, y1));
  }
  for (int y = y1; y > y0; --y) {
    p.push_back(Point(x0, y));
  }
  for (int x = x0; x < x_mid; ++x) {
    p.push_back(Point(x, y0));
  }
  return p;
}

struct Division {
  std::string id;
  std::vector<Polygon_with_holes> pwhs;

  // Ratio of target area to area; distinct ratios give the divisions
  // distinct target densities
  double target_ratio;
};

// Exact overlap area of a ring and the cell [x, x + 1] x [y, y + 1]. We clip
// the ring by each side of the cell with the Sutherland-Hodgman algorithm,
// which gives the exact overlap area because the cell is convex.
double overlap_area(
  const Polygon &r,
  const unsigned int x,
  const unsigned int y)
{
  std::vector<Point> pts(r.vertices_begin(), r.vertices_end());

  // Keep the points with sign * (coordinate - value) >= 0, where the
  // coordinate is x if `along_x` is true and y otherwise
  const auto clip =
    [&pts](const bool along_x, const double value, const double sign) {
      const auto dist = [&](const Point &p) {
        return sign * ((along_x ? p.x() : p.y()) - value);
      };
      std::vector<Point> out;
      for (std::size_t i = 0; i < pts.size(); ++i) {
        const Point &p = pts[(i + pts.size() - 1) % pts.size()];
        const Point &q = pts[i];
        const double dp = dist(p);
        const double dq = dist(q);
        if ((dp >= 0) != (dq >= 0)) {
          const double t = dp / (dp - dq);
          out.emplace_back(
            p.x() + t * (q.x() - p.x()),
            p.y() + t * (q.y() - p.y()));
        }
        if (dq >= 0) {
          out.push_back(q);
        }
      }
      pts = out;
    };
  clip(true, x, 1.0);
  clip(true, x + 1.0, -1.0);
  clip(false, y, 1.0);
  clip(false, y + 1.0, -1.0);
  double twice_area = 0.0;
  for (std::size_t i = 0; i < pts.size(); ++i) {
    const Point &p = pts[i];
    const Point &q = pts[(i + 1) % pts.size()];
    twice_area += p.x() * q.y() - q.x() * p.y();
  }
  return std::abs(twice_area) / 2.0;
}

double overlap_area(
  const Polygon_with_holes &pwh,
  const unsigned int x,
  const unsigned int y)
{
  double a = overlap_area(pwh.outer_boundary(), x, y);
  for (const auto &h : pwh.holes()) {
    a -= overlap_area(h, x, y);
  }
  return a;
}

// Fill the lattice with the density of the divisions and compare each cell
// with the expected density. Return the sum of the densities of all cells.
double check_density(const std::vector<Division> &divisions)
{
  Arguments args{};
  InsetState inset("C", args);
  for (const auto &d : divisions) {
    GeoDiv gd(d.id);
    for (const auto &pwh : d.pwhs) {
      gd.push_back(pwh);
    }
    inset.push_back(gd);
    inset.insert_target_area(d.id, d.target_ratio * gd.area());
  }
  inset.set_grid_dimensions(lx, ly);
  inset.store_initial_area();
  inset.set_area_errors();
  inset.ref_to_rho_init().allocate(lx, ly);
  inset.ref_to_rho_ft().allocate(lx, ly);
  inset.make_fftw_plans_for_rho();
  inset.fill_with_density();

  // Weighted average as documented above
  const double n_cells = lx * ly;
  const double inset_area = inset.total_inset_area();
  const double target_area = inset.total_target_area();
  const double ocean_density =
    (n_cells - target_area) / (n_cells - inset_area);
  const double ocean_area_error =
    std::abs((n_cells - inset_area) / (n_cells - target_area) - 1.0);
  std::ostringstream wrong_cells;
  unsigned int n_wrong = 0;
  double max_diff = 0.0;
  double sum_rho = 0.0;
  for (unsigned int i = 0; i < lx; ++i) {
    for (unsigned int j = 0; j < ly; ++j) {
      double num = 0.0;
      double den = 0.0;
      double covered = 0.0;
      for (std::size_t k = 0; k < divisions.size(); ++k) {
        const GeoDiv &gd = inset.geo_divs()[k];
        double a = 0.0;
        for (const auto &pwh : gd.polygons_with_holes()) {
          a += overlap_area(pwh, i, j);
        }
        const double w = a * inset.area_error_at(gd.id());
        num += w * inset.target_area_at(gd.id()) / gd.area();
        den += w;
        covered += a;
      }
      const double w_ocean = (1.0 - covered) * ocean_area_error;
      num += w_ocean * ocean_density;
      den += w_ocean;
      const double expected = den > 0.0 ? num / den : ocean_density;
      const double rho = inset.ref_to_rho_init()(i, j);
      sum_rho += rho;
      const double diff = std::abs(rho - expected);
      max_diff = std::max(max_diff, diff);
      if (diff > 1e-12) {
        ++n_wrong;
        wrong_cells << " (" << i << ", " << j << "): " << rho << " instead of "
                    << expected << ";";
      }
    }
  }
  BOOST_TEST_MESSAGE("max |rho - expected| = " << max_diff);
  BOOST_TEST(
    n_wrong == 0u,
    n_wrong << " cells with wrong density:" << wrong_cells.str());
  inset.destroy_fftw_plans_for_rho();
  inset.ref_to_rho_init().free();
  inset.ref_to_rho_ft().free();
  return sum_rho;
}

}  // namespace

// The rectangles cover whole cells, so that the density of each cell is
// either the target density of a rectangle or the ocean density. The sum of
// the densities is then the total target area plus the ocean area times the
// ocean density, which is lx * ly.
BOOST_AUTO_TEST_CASE(rectangles_with_sides_on_lattice_lines)
{
  const std::vector<Division> divisions{
    {"plain",
     {Polygon_with_holes(ring({{4, 3}, {12, 3}, {12, 9}, {4, 9}}))},
     2.0},
    {"densified",
     {Polygon_with_holes(densified_rectangle(17, 20, 27, 28))},
     0.5}};
  const double sum_rho = check_density(divisions);
  BOOST_TEST(std::abs(sum_rho - lx * ly) < 1e-9);
}

// Diagonal edges through lattice points cut cells along their diagonals and
// touch other cells at a corner only. The notch of the second polygon has a
// reflex vertex at (6, 24), which touches the corners of cells inside the
// polygon.
BOOST_AUTO_TEST_CASE(diagonal_edges_through_lattice_points)
{
  const std::vector<Division> divisions{
    {"diamond",
     {Polygon_with_holes(ring({{24, 4}, {28, 8}, {24, 12}, {20, 8}}))},
     2.0},
    {"notch",
     {Polygon_with_holes(
       ring({{2, 20}, {10, 20}, {10, 28}, {6, 24}, {2, 28}}))},
     0.5}};
  check_density(divisions);
}

// Holes are clockwise. The cells inside the holes must get the ocean density.
BOOST_AUTO_TEST_CASE(holes_with_edges_on_lattice_lines)
{
  Polygon_with_holes square_hole(ring({{2, 2}, {14, 2}, {14, 14}, {2, 14}}));
  square_hole.add_hole(ring({{5, 6}, {5, 10}, {9, 10}, {9, 6}}));
  Polygon_with_holes diamond_hole(
    ring({{16, 16}, {30, 16}, {30, 30}, {16, 30}}));
  diamond_hole.add_hole(ring({{23, 18}, {19, 22}, {23, 26}, {27, 22}}));
  const std::vector<Division> divisions{
    {"square hole", {square_hole}, 2.0},
    {"diamond hole", {diamond_hole}, 0.5}};
  check_density(divisions);
}

// Neighboring divisions share edges on lattice lines and a vertex at the
// lattice point (10, 10). Each cell must be counted once, so the sum of the
// densities is lx * ly, as for the rectangles above.
BOOST_AUTO_TEST_CASE(divisions_sharing_edges_on_lattice_lines)
{
  const std::vector<Division> divisions{
    {"A",
     {Polygon_with_holes(ring({{4, 4}, {10, 4}, {10, 10}, {4, 10}}))},
     2.0},
    {"B",
     {Polygon_with_holes(ring({{10, 4}, {16, 4}, {16, 10}, {10, 10}}))},
     0.5},
    {"C",
     {Polygon_with_holes(
       ring({{4, 10}, {10, 10}, {16, 10}, {16, 14}, {4, 14}}))},
     3.0}};
  const double sum_rho = check_density(divisions);
  BOOST_TEST(std::abs(sum_rho - lx * ly) < 1e-9);
}

// Two triangles share a diagonal through lattice points. The cells along the
// diagonal are split in half between the triangles.
BOOST_AUTO_TEST_CASE(divisions_sharing_a_diagonal_through_lattice_points)
{
  const std::vector<Division> divisions{
    {"below", {Polygon_with_holes(ring({{20, 2}, {28, 2}, {28, 10}}))}, 2.0},
    {"above", {Polygon_with_holes(ring({{20, 2}, {28, 10}, {20, 10}}))}, 0.5}};
  check_density(divisions);
}
