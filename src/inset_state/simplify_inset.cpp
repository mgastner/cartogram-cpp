// TODO: What happens if two polygons have touching lines but the corner
// points are not identical in both lines?

#include "../inset_state.h"
#include "../constants.h"
#include <algorithm>

// We use -1 to signal that there is no simplified polygon that can be matched
// with a given non-simplified polygon
constexpr int no_matching_simpl_pgn = -1;

bool contains_vertices_in_order(const Polygon non_simpl_pgn,
                                const Polygon simpl_pgn,
                                const Bbox simpl_bb)
{
  // Return true if and only if
  // - the non-simplified polygon contains all vertices in the simplified
  //   polygon and
  // - the order of the vertices in both polygons match.
  // We pass the bounding box of the simplified polygon as an argument so that
  // we do not need to construct the bounding box in this function.
  if (non_simpl_pgn.size() < simpl_pgn.size()) {
    return false;  // Simplification cannot create additional vertices
  }

  // If the two bounding boxes neither touch nor overlap, then the
  // non-simplified polygon does not contain the vertices of the simplified
  // polygon. For the if-condition, see:
  // https://stackoverflow.com/questions/325933/determine-whether-two-date-
  // ranges-overlap/325964#325964 (accessed on 2021-10-05).
  const auto non_simpl_bb = non_simpl_pgn.bbox();
  if (non_simpl_bb.xmax() < simpl_bb.xmin() ||
      non_simpl_bb.xmin() > simpl_bb.xmax() ||
      non_simpl_bb.ymax() < simpl_bb.ymin() ||
      non_simpl_bb.ymin() > simpl_bb.ymax()) {
    return false;
  }
  std::vector<unsigned int> indices;
  for (const auto &simpl_pt : simpl_pgn) {
    const auto non_simpl_it = std::find(non_simpl_pgn.vertices_begin(),
                                        non_simpl_pgn.vertices_end(),
                                        simpl_pt);

    // Return false if there is no matching vertex in the non-simplified
    // polygon
    if (non_simpl_it == non_simpl_pgn.vertices_end()) {
      return false;
    }
    indices.push_back(distance(non_simpl_pgn.vertices_begin(), non_simpl_it));
    if (!std::is_sorted(indices.begin(), indices.end())) {
      return false;
    }
  }
  return std::is_sorted(indices.begin(), indices.end());
}

int simplified_polygon_index(const Polygon non_simpl_pgn,
                             const std::vector<Polygon> *simpl_pgns,
                             const std::vector<Bbox> *simpl_bboxes,
                             std::list<unsigned int> *unmatched)
{
  // Return index of simplified polygon that corresponds to a non-simplified
  // polygon. We pass a vector of the bounding boxes of the simplified
  // polygons as an argument so that the bounding boxes do not need to be
  // calculated repeatedly. We maintain a list of hitherto unmatched polygon
  // indices to avoid unnecessary iterations.
  for (const auto &i : *unmatched) {
    if (contains_vertices_in_order(non_simpl_pgn,
                                   simpl_pgns->at(i),
                                   simpl_bboxes->at(i))) {
      unmatched->remove(i);
      return i;
    }
  }
  return no_matching_simpl_pgn;
}

void InsetState::simplify(const unsigned int target_points_per_inset)
{
  const unsigned int n_pts_before = n_points();
  std::cerr << n_pts_before << " points in inset. ";

  if (n_pts_before <= target_points_per_inset) {
    std::cerr << "No need for simplification." << std::endl;
    return;
  }
  std::cerr << "Simplifying the inset. " << std::endl;

  // Store Polygons as a CT (Constrained Triangulation) object. Code inspired
  // by https://doc.cgal.org/latest/Polyline_simplification_2/index.html
  CT ct;
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      ct.insert_constraint(pwh.outer_boundary());
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        ct.insert_constraint(*h);
      }
    }
  }

  // Simplify polygons
  const unsigned long target_pts =
    std::max(target_points_per_inset,
             min_points_per_ring * n_rings());
  const double ratio = static_cast<double>(target_pts) / n_pts_before;
  PS::simplify(ct, Cost(), Stop(ratio));

  // Store each constraint in ct as a polygon. Also store bounding box so
  // that we can match non-simplified and simplified polygons more quickly.
  std::vector<Polygon> simpl_pgns;
  std::vector<Bbox> simpl_bboxes;
  for (auto it = ct.constraints_begin(); it != ct.constraints_end(); ++it) {

    // First and last point in the constraint are identical. We remove the
    // last point to make the polygon simple.
    const Polygon ct_as_pgn(ct.points_in_constraint_begin(*it),
                            --ct.points_in_constraint_end(*it));
    simpl_pgns.push_back(ct_as_pgn);
    simpl_bboxes.push_back(ct_as_pgn.bbox());
  }

  // Match non-simplified polygon to its simplified counterpart
  std::list<unsigned int> unmatched(
    boost::counting_iterator<unsigned int>(0U),
    boost::counting_iterator<unsigned int>(simpl_pgns.size()));
  std::vector<int> matching_simpl_pgn;
  for (const auto &gd : geo_divs_) {
    for (const auto &pwh : gd.polygons_with_holes()) {
      const int ext_index = simplified_polygon_index(pwh.outer_boundary(),
                                                     &simpl_pgns,
                                                     &simpl_bboxes,
                                                     &unmatched);
      matching_simpl_pgn.push_back(ext_index);
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        const int hole_index = simplified_polygon_index(*h,
                                                        &simpl_pgns,
                                                        &simpl_bboxes,
                                                        &unmatched);
        matching_simpl_pgn.push_back(hole_index);
      }
    }
  }

  // Sanity check
  if (std::find(matching_simpl_pgn.begin(),
                matching_simpl_pgn.end(),
                no_matching_simpl_pgn) != matching_simpl_pgn.end() ||
      !unmatched.empty()) {
    std::cerr << "ERROR: Unmatched polygon in "
              << __func__
              << "()."
              << std::endl;
    for (const auto &u : unmatched) {
      std::cerr << "Unmatched polygon: "
                << u
                << std::endl;
      for (const auto &pt : simpl_pgns[u]) {
        std::cerr << pt << std::endl;
      }
    }
    exit(1);
  }

  // Replace non-simplified polygons by their simplified counterparts
  unsigned int pgn_ctr = 0;
  for (auto &gd : geo_divs_) {
    for (auto &pwh : *gd.ref_to_polygons_with_holes()) {
      const unsigned int match = matching_simpl_pgn[pgn_ctr++];
      pwh.outer_boundary() = simpl_pgns[match];
      for (auto h = pwh.holes_begin(); h != pwh.holes_end(); ++h) {
        const unsigned int match = matching_simpl_pgn[pgn_ctr++];
        *h = simpl_pgns[match];
      }
    }
  }
  std::cerr << n_points() << " points after simplification." << std::endl;
  return;
}
