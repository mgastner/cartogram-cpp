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
  std::cerr << "Simplifying the inset. " << std::endl;

  // Store Polygons as a CT (Constrained Triangulation) object. Code inspired
  // by https://doc.cgal.org/latest/Polyline_simplification_2/index.html
  std::unordered_map<int, Constraint_id> pgn_id_to_constraint_id;
  int pgn_id = 0;
  CT ct;
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      pgn_id_to_constraint_id[pgn_id++] =
        // NOTE: This line causes a segmentation fault in
        // cartogram ../../sample_data/world_by_country_since_2022/world_by_country_since_2022.geojson ../../sample_data/world_by_country_since_2022/world_population_2010.csv -QST
        // in DEBUG mode
        // TODO: Investigate why
        ct.insert_constraint(pwh.outer_boundary());
      for (const auto &h : pwh.holes()) {
        pgn_id_to_constraint_id[pgn_id++] = ct.insert_constraint(h);
      }
    }
  }

  // Simplify polygons
  const unsigned long target_pts =
    std::max(target_points_per_inset, min_points_per_ring * n_rings());
  const double ratio = static_cast<double>(target_pts) / n_pts_before;
  CGAL::Polyline_simplification_2::simplify(ct, Cost(), Stop(ratio));

  // Store simplified polygons
  pgn_id = 0;
  for (auto &gd : geo_divs_) {
    for (auto &pwh : gd.ref_to_polygons_with_holes()) {
      auto cit = pgn_id_to_constraint_id[pgn_id++];

      Polygon ext_ring;
      for (auto it = ct.vertices_in_constraint_begin(cit);
           it != ct.vertices_in_constraint_end(cit);
           ++it) {
        ext_ring.push_back(((*it)->point()));
      }

      // remove the last point which is identical to the first point to
      // make it a simple polygon
      ext_ring.erase(ext_ring.vertices_end() - 1);

      std::vector<Polygon> int_ring_v;
      for (auto h : pwh.holes()) {
        cit = pgn_id_to_constraint_id[pgn_id++];
        Polygon int_ring;
        for (auto it = ct.vertices_in_constraint_begin(cit);
             it != ct.vertices_in_constraint_end(cit);
             ++it) {
          int_ring.push_back(((*it)->point()));
        }

        // remove the last point which is identical to the first point to
        // make it a simple polygon
        int_ring.erase(int_ring.vertices_end() - 1);
        int_ring_v.push_back(int_ring);
      }
      Polygon_with_holes pgnwh =
        Polygon_with_holes(ext_ring, int_ring_v.begin(), int_ring_v.end());
      pwh = std::move(pgnwh);
    }
  }

  std::cerr << n_points() << " points after simplification." << std::endl;
  timer.stop("Simplification");
}
