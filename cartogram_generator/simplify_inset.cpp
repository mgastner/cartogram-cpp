// TODO: What happens if two polygons have touching lines but the corner
// points are not identical in both lines?

// TODO: typedef CGAL::BBox_2 as Bbox

#include "constants.h"
#include "inset_state.h"
#include <algorithm>

// We use the value no_matching_non_simpl_pgon to signal that there is no
// simplified polygon that we could match with a non-simplified polygon
constexpr int no_matching_non_simpl_pgon = -1;

bool contains_vertices_in_order(const Polygon non_simpl_pgon,
                                const Polygon simpl_pgon,
                                const CGAL::Bbox_2 simpl_bbox)
{
  // Return true if and only if the non-simplified polygon contains all
  // vertices in the simplified polygon and the order of the vertices in
  // both polygons match. We pass the bounding box of the simplified bounding
  // box as a parameter so that we do not need to construct it in the function
  // body.
  if (non_simpl_pgon.size() < simpl_pgon.size()) {
    return false;  // Simplification cannot create additional vertices
  }

  // If the two bounding boxes neither touch nor overlap, then the
  // non-simplified polygon does not contain the vertices of the simplified
  // polygon. For the if-condition, see:
  // https://stackoverflow.com/questions/325933/determine-whether-two-date-
  // ranges-overlap/325964#325964 (accessed on 2021-10-05).
  CGAL::Bbox_2 non_simpl_bbox = non_simpl_pgon.bbox();
  if (non_simpl_bbox.xmax() < simpl_bbox.xmin()
      || non_simpl_bbox.xmin() > simpl_bbox.xmax()
      || non_simpl_bbox.ymax() < simpl_bbox.ymin()
      || non_simpl_bbox.ymin() > simpl_bbox.ymax()) {
    return false;
  }
  std::vector<unsigned int> indices;
  for (const auto simpl_pt : simpl_pgon) {
    const auto non_simpl_it = std::find(non_simpl_pgon.vertices_begin(),
                                        non_simpl_pgon.vertices_end(),
                                        simpl_pt);

    // Return false if there is no matching vertex in polygon1
    if (non_simpl_it == non_simpl_pgon.vertices_end()) {
      return false;
    }
    indices.push_back(distance(non_simpl_pgon.vertices_begin(),
                               non_simpl_it));
    if (!std::is_sorted(indices.begin(), indices.end())) {
      return false;
    }
  }
  return std::is_sorted(indices.begin(), indices.end());
}

int simplified_polygon_index(const Polygon non_simpl_pgon,
                             const std::vector<Polygon> *simpl_pgons,
                             const std::vector<CGAL::Bbox_2> *simpl_bboxes,
                             std::list<unsigned int> *unmatched)
{
  // Return index of polygon in simpl_pgon that corresponds to a
  // non-simplified polygon. We pass the bounding boxes of the simplified
  // polygons as a parameter so that the bounding boxes do not need to be
  // calculated when calling contains_vertices_in_order(). We maintain a list
  // of hitherto unmatched polygon indices to avoid unnecessary calls of
  // contains_vertices_in_order().
  for (const auto i : *unmatched) {
    if (contains_vertices_in_order(non_simpl_pgon,
                                   simpl_pgons->at(i),
                                   simpl_bboxes->at(i))) {
      unmatched->remove(i);
      return i;
    }
  }
  return no_matching_non_simpl_pgon;
}

void simplify_inset(InsetState *inset_state)
{
  const unsigned int n_pts_before = inset_state->n_points();
  if (n_pts_before <= target_points_per_inset) {
    std::cerr << n_pts_before
              << " points in inset. No need for simplification."
              << std::endl;
    return;
  }
  std::cerr << "Simplifying the inset. "
            << n_pts_before
            << " points in the inset before simplification."
            << std::endl;

  // Store Polygons as CT (Constrained Triangulation) objects. Code inspired
  // by https://doc.cgal.org/latest/Polyline_simplification_2/index.html
  CT ct;
  for (const auto gd : inset_state->geo_divs()) {
    for (const auto pwh : gd.polygons_with_holes()) {
      ct.insert_constraint(pwh.outer_boundary());
      for (auto hi = pwh.holes_begin(); hi != pwh.holes_end(); ++hi) {
        ct.insert_constraint(*hi);
      }
    }
  }

  // Simplify polygons
  const double ratio = double(target_points_per_inset) / n_pts_before;
  PS::simplify(ct, Cost(), Stop(ratio));

  // Store each constraint in ct as a polygon. Also store bounding box so
  // that we can match non-simplified and simplified polygons more quickly.
  std::vector<Polygon> simpl_pgons;
  std::vector<CGAL::Bbox_2> simpl_bboxes;
  for (auto it = ct.constraints_begin(); it != ct.constraints_end(); ++it) {

    // First and last point in the constraint are identical. We remove the
    // last point to make the polygon simple.
    const Polygon ct_as_pgon(ct.points_in_constraint_begin(*it),
                             --ct.points_in_constraint_end(*it));
    simpl_pgons.push_back(ct_as_pgon);
    simpl_bboxes.push_back(ct_as_pgon.bbox());
  }

  // Match non-simplified polygon to its simplified counterpart
  std::list<unsigned int> unmatched(
    boost::counting_iterator<unsigned int>(0U),
    boost::counting_iterator<unsigned int>(simpl_pgons.size()));
  std::vector<int> matching_simpl_pgon;
  for (const auto gd : inset_state->geo_divs()) {
    for (const auto pwh : gd.polygons_with_holes()) {
      const int ext_index = simplified_polygon_index(pwh.outer_boundary(),
                                                     &simpl_pgons,
                                                     &simpl_bboxes,
                                                     &unmatched);
      matching_simpl_pgon.push_back(ext_index);
      for (auto hi = pwh.holes_begin(); hi != pwh.holes_end(); ++hi) {
        const int hole_index = simplified_polygon_index(*hi,
                                                        &simpl_pgons,
                                                        &simpl_bboxes,
                                                        &unmatched);
        matching_simpl_pgon.push_back(hole_index);
      }
    }
  }

  // Sanity check
  if (std::find(matching_simpl_pgon.begin(),
                matching_simpl_pgon.end(),
                no_matching_non_simpl_pgon) != matching_simpl_pgon.end() ||
      !unmatched.empty()) {
    std::cerr << "Unmatched polygon in simplify_inset()." << std::endl;
    exit(1);
  }

  // Replace non-simplified polygons by their simplified counterparts
  unsigned int pgon_ctr = 0;
  for (auto &gd : *inset_state->ref_to_geo_divs()) {
    for (auto &pwh : *gd.ref_to_polygons_with_holes()) {
      const unsigned int match = matching_simpl_pgon[pgon_ctr++];
      pwh.outer_boundary() = simpl_pgons[match];
      for (auto hi = pwh.holes_begin(); hi != pwh.holes_end(); ++hi) {
        const unsigned int match = matching_simpl_pgon[pgon_ctr++];
        *hi = simpl_pgons[match];
      }
    }
  }
  std::cerr << inset_state->n_points()
            << " points after simplification."
            << std::endl;
  return;
}
