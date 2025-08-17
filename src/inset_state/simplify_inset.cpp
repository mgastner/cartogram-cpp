// TODO: What happens if two polygons have touching lines, but the corner
//       points are not identical in both lines?

#include "constants.hpp"
#include "inset_state.hpp"

using Vb = CGAL::Polyline_simplification_2::Vertex_base_2<Scd>;
using Fb = CGAL::Constrained_triangulation_face_base_2<Scd>;
using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using CDT = CGAL::
  Constrained_Delaunay_triangulation_2<Scd, TDS, CGAL::Exact_predicates_tag>;
using CT = CGAL::Constrained_triangulation_plus_2<CDT>;
using Constraint_id = CT::Constraint_id;
using Stop = CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold;
using Cost = CGAL::Polyline_simplification_2::Squared_distance_cost;

void InsetState::simplify(const unsigned int target_points_per_inset)
{
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
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      pgn_id_to_constraint_id[pgn_id++] =
        ct.insert_constraint(pwh.outer_boundary());
      for (const auto &h : pwh.holes()) {
        pgn_id_to_constraint_id[pgn_id++] = ct.insert_constraint(h);
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
        ext_ring.push_back(*pit);

      std::vector<Polygon> int_ring_v;
      int_ring_v.reserve(pwh.number_of_holes());
      for (auto h : pwh.holes()) {
        cit = pgn_id_to_constraint_id[pgn_id++];
        Polygon int_ring;
        for (auto pit = ct.points_in_constraint_begin(cit),
                  pend = ct.points_in_constraint_end(cit);
             std::next(pit) != pend;
             ++pit)
          int_ring.push_back(*pit);
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
