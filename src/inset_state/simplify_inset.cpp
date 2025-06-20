// TODO: What happens if two polygons have touching lines, but the corner
//       points are not identical in both lines?

#include "constants.hpp"
#include "inset_state.hpp"

void InsetState::simplify(const unsigned int target_points_per_inset)
{
  timer.start("Simplification");
  const size_t n_pts_before = n_points();
  std::cerr << n_pts_before << " points in inset. ";
  if (n_pts_before <= target_points_per_inset) {
    std::cerr << "No need for simplification." << std::endl;
    return;
  }
  const unsigned int n_rings = this->n_rings();
  const unsigned long target_pts =
    std::max(target_points_per_inset, min_points_per_ring * n_rings);
  const double ratio =
    static_cast<double>(target_pts) / static_cast<double>(n_pts_before);
  std::cerr << "Target points: " << target_pts << std::endl;
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
      auto beg = ct.vertices_in_constraint_begin(cit);
      auto end = ct.vertices_in_constraint_end(cit);

      Polygon ext_ring;
      const std::ptrdiff_t len = std::distance(beg, end);
      if (len > 1)
        ext_ring.reserve(static_cast<std::size_t>(len - 1));
      for (auto it = beg; std::next(it) != end; ++it)
        ext_ring.push_back((*it)->point());

      std::vector<Polygon> int_ring_v;
      int_ring_v.reserve(pwh.number_of_holes());
      for (auto h : pwh.holes()) {
        cit = pgn_id_to_constraint_id[pgn_id++];
        beg = ct.vertices_in_constraint_begin(cit);
        end = ct.vertices_in_constraint_end(cit);

        Polygon int_ring;
        const std::ptrdiff_t ring_len = std::distance(beg, end);
        if (ring_len > 1)
          ext_ring.reserve(static_cast<std::size_t>(ring_len - 1));
        for (auto it = beg; std::next(it) != end; ++it)
          int_ring.push_back((*it)->point());

        int_ring_v.emplace_back(std::move(int_ring));
      }

      Polygon_with_holes pgnwh(
        std::move(ext_ring),
        int_ring_v.begin(),
        int_ring_v.end());
      pwh = std::move(pgnwh);
    }
  }

  std::cerr << n_points() << " points after simplification." << std::endl;
  timer.stop("Simplification");
}
