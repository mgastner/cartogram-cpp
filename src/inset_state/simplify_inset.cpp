#include "constants.hpp"
#include "inset_state.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// The constrained triangulation decides whether vertices of one ring lie on
// another ring's segments and whether segments cross. These decisions need
// exact predicates: with plain double arithmetic, nearly collinear vertices
// on a shared border are misclassified, the copies of the border are merged
// only in part, and spurious intersection points are constructed. The
// coordinates themselves remain doubles. Borders that two GeoDivs share are
// simplified identically if both copies have the same vertices. The repair
// of the input topology (topology_repair.hpp) gives both copies the same
// vertices where they overlap or miss each other by less than the snapping
// distance, unless --disable_topology_repair is passed.
using SimplifyKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Polyline_simplification_2::Vertex_base_2<SimplifyKernel>;
using Fb = CGAL::Constrained_triangulation_face_base_2<SimplifyKernel>;
using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<
  SimplifyKernel,
  TDS,
  CGAL::Exact_predicates_tag>;
using CT = CGAL::Constrained_triangulation_plus_2<CDT>;
using Constraint_id = CT::Constraint_id;
using Stop = CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold;
using Cost = CGAL::Polyline_simplification_2::Squared_distance_cost;

void InsetState::simplify(const unsigned int target_points_per_inset)
{
  if (intersections_found_) {
    std::cerr
      << "WARNING: Intersections found previously; skipping simplification!"
      << std::endl;
    return;
  }

  timer.start("Simplification");
  const size_t n_pts_before = n_points();
  std::cerr << n_pts_before << " points in inset. ";

  const unsigned int n_rings = this->n_rings();
  const unsigned long target_pts =
    std::max(target_points_per_inset, min_points_per_ring * n_rings);
  std::cerr << "Target points: " << target_pts << std::endl;

  if (n_pts_before <= target_pts) {
    std::cerr << "No need for simplification." << std::endl;
    return;
  }
  const double ratio =
    static_cast<double>(target_pts) / static_cast<double>(n_pts_before);
  std::cerr << "Simplifying the inset. " << std::endl;

  // Store Polygons as a CT (Constrained Triangulation) object. Code inspired
  // by https://doc.cgal.org/latest/Polyline_simplification_2/index.html
  std::vector<Constraint_id> pgn_id_to_constraint_id(n_rings);
  size_t pgn_id = 0;
  CT ct;
  auto insert_ring = [&ct](const Polygon &ring) {
    std::vector<CT::Point> pts;
    pts.reserve(ring.size());
    for (const auto &p : ring) {
      pts.emplace_back(p.x(), p.y());
    }
    return ct.insert_constraint(
      CGAL::Polygon_2<SimplifyKernel>(pts.begin(), pts.end()));
  };
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      pgn_id_to_constraint_id[pgn_id++] = insert_ring(pwh.outer_boundary());
      for (const auto &h : pwh.holes()) {
        pgn_id_to_constraint_id[pgn_id++] = insert_ring(h);
      }
    }
  }

  // Simplify polygons
  CGAL::Polyline_simplification_2::simplify(ct, Cost(), Stop(ratio));

  // Store simplified polygons
  pgn_id = 0;
  for (auto &gd : geo_divs_) {
    for (auto &pwh : gd.ref_to_polygons_with_holes()) {
      auto cit = pgn_id_to_constraint_id[pgn_id++];

      Polygon ext_ring;
      for (auto pit = ct.points_in_constraint_begin(cit),
                pend = ct.points_in_constraint_end(cit);
           std::next(pit) != pend;
           ++pit)
        ext_ring.push_back(Point(pit->x(), pit->y()));

      std::vector<Polygon> int_ring_v;
      int_ring_v.reserve(pwh.number_of_holes());
      for (auto h : pwh.holes()) {
        cit = pgn_id_to_constraint_id[pgn_id++];
        Polygon int_ring;
        for (auto pit = ct.points_in_constraint_begin(cit),
                  pend = ct.points_in_constraint_end(cit);
             std::next(pit) != pend;
             ++pit)
          int_ring.push_back(Point(pit->x(), pit->y()));
        int_ring_v.emplace_back(std::move(int_ring));
      }

      pwh = Polygon_with_holes(
        std::move(ext_ring),
        int_ring_v.begin(),
        int_ring_v.end());
    }
  }

  std::cerr << n_points() << " points after simplification." << std::endl;
  timer.stop("Simplification");
}
