#include "inset_state.hpp"
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <cmath>

// With --periodic_x, regions flow across the seam at x = 0 and x = lx, so
// that a polygon can extend past the sides of the frame. Before the output,
// each such polygon is cut at the lines x = k * lx, and each piece is shifted
// by -k * lx into the frame. The cuts use exact constructions, so that the
// pieces on the two sides of a cut, and neighbors that share a border across
// it, get the same cut points. The exact Boolean operations of CGAL require
// valid polygons, so a polygon that is not valid is clipped ring by ring
// instead. The unsimplified input geometry that --redirect_exports_to_stdout
// writes along is not cut: go-cart.io morphs its vertices one to one into
// the vertices of the input map, so its coordinates may extend past the
// frame.

using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;
using Exact_polygon = CGAL::Polygon_2<Epeck>;
using Exact_polygon_with_holes = CGAL::Polygon_with_holes_2<Epeck>;

static Exact_polygon exact_ring(const Polygon &ring)
{
  Exact_polygon exact;
  for (const Point &p : ring) {
    exact.push_back(Epeck::Point_2(p.x(), p.y()));
  }
  return exact;
}

// The ring shifted by -shift in x, rounded to doubles
static Polygon rounded_ring(const Exact_polygon &ring, const double shift)
{
  Polygon rounded;
  for (const auto &p : ring) {
    rounded.push_back(Point(
      CGAL::to_double(CGAL::exact(p.x() - shift)),
      CGAL::to_double(CGAL::exact(p.y()))));
  }
  return rounded;
}

// The part of the ring in the strip x_min <= x <= x_max, shifted by -x_min
// (Sutherland-Hodgman clipping against the two sides of the strip). Where the
// ring leaves and reenters the strip, the result runs along its side.
static Polygon clipped_ring(
  const Polygon &ring,
  const double x_min,
  const double x_max)
{
  std::vector<Point> pts(ring.begin(), ring.end());
  for (const int side : {0, 1}) {
    const double c = (side == 0) ? x_min : x_max;
    auto inside = [&](const Point &p) {
      return (side == 0) ? !(p.x() < c) : !(p.x() > c);
    };
    std::vector<Point> out;
    for (std::size_t k = 0; k < pts.size(); ++k) {
      const Point &a = pts[k];
      const Point &b = pts[(k + 1) % pts.size()];
      if (inside(a)) {
        out.push_back(a);
      }
      if (inside(a) != inside(b)) {
        const double t = (c - a.x()) / (b.x() - a.x());
        out.emplace_back(c, a.y() + t * (b.y() - a.y()));
      }
    }
    pts = std::move(out);
    if (pts.empty()) {
      break;
    }
  }
  Polygon clipped;
  for (const Point &p : pts) {
    clipped.push_back(Point(p.x() - x_min, p.y()));
  }
  return clipped;
}

// Appends the pieces of pwh, each shifted into the frame 0 <= x <= lx, to
// pieces. Returns false if pwh is not valid and was clipped ring by ring.
static bool wrap_polygon_with_holes(
  const Polygon_with_holes &pwh,
  const double lx,
  std::vector<Polygon_with_holes> &pieces)
{
  const auto bb = pwh.outer_boundary().bbox();
  if (!(bb.xmin() < 0.0) && !(bb.xmax() > lx)) {
    pieces.push_back(pwh);
    return true;
  }
  Exact_polygon outer = exact_ring(pwh.outer_boundary());
  if (outer.is_clockwise_oriented()) {
    outer.reverse_orientation();
  }
  std::vector<Exact_polygon> holes;
  for (const auto &h : pwh.holes()) {
    holes.push_back(exact_ring(h));
    if (holes.back().is_counterclockwise_oriented()) {
      holes.back().reverse_orientation();
    }
  }
  const Exact_polygon_with_holes exact(outer, holes.begin(), holes.end());

  // Strips k * lx <= x <= (k + 1) * lx that the polygon reaches
  const auto k_min = static_cast<long>(std::floor(bb.xmin() / lx));
  const auto k_max = static_cast<long>(std::ceil(bb.xmax() / lx)) - 1;
  const double y_min = bb.ymin() - 1.0;
  const double y_max = bb.ymax() + 1.0;
  const bool valid = CGAL::is_valid_polygon_with_holes(
    exact,
    CGAL::Gps_segment_traits_2<Epeck>());
  for (long k = k_min; k <= k_max; ++k) {
    const double x_min = static_cast<double>(k) * lx;
    const double x_max = static_cast<double>(k + 1) * lx;
    if (!valid) {
      const Polygon piece = clipped_ring(pwh.outer_boundary(), x_min, x_max);
      if (piece.size() < 3) {
        continue;
      }
      std::vector<Polygon> piece_holes;
      for (const auto &h : pwh.holes()) {
        const Polygon clipped = clipped_ring(h, x_min, x_max);
        if (clipped.size() >= 3) {
          piece_holes.push_back(clipped);
        }
      }
      pieces.emplace_back(piece, piece_holes.begin(), piece_holes.end());
      continue;
    }
    Exact_polygon strip;
    strip.push_back(Epeck::Point_2(x_min, y_min));
    strip.push_back(Epeck::Point_2(x_max, y_min));
    strip.push_back(Epeck::Point_2(x_max, y_max));
    strip.push_back(Epeck::Point_2(x_min, y_max));
    std::vector<Exact_polygon_with_holes> parts;
    CGAL::intersection(exact, strip, std::back_inserter(parts));
    for (const auto &part : parts) {
      std::vector<Polygon> part_holes;
      for (const auto &h : part.holes()) {
        part_holes.push_back(rounded_ring(h, x_min));
      }
      pieces.emplace_back(
        rounded_ring(part.outer_boundary(), x_min),
        part_holes.begin(),
        part_holes.end());
    }
  }
  return valid;
}

void InsetState::wrap_into_frame()
{
  if (!args_.periodic_x) {
    return;
  }
  const double dlx = lx_;
  unsigned int n_invalid = 0;
  for (auto &gd : geo_divs_) {
    std::vector<Polygon_with_holes> pieces;
    for (const auto &pwh : gd.polygons_with_holes()) {
      if (!wrap_polygon_with_holes(pwh, dlx, pieces)) {
        ++n_invalid;
      }
    }
    gd.ref_to_polygons_with_holes() = std::move(pieces);
  }
  if (n_invalid > 0) {
    std::cerr << "WARNING: " << n_invalid << " polygon(s) across the seam "
              << "are not valid; they were clipped ring by ring." << std::endl;
  }
}
