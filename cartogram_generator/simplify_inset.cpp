// TODO: What happens if two polygons have touching lines but the corner
// points are not identical in both lines?

// TODO: typedef CGAL::BBox_2 as Bbox

#include <algorithm>
#include "inset_state.h"

constexpr int no_matching_non_simplified_polygon = -1;

bool contains_vertices_in_order(Polygon non_simplified_polygon,
                                Polygon simplified_polygon,
                                CGAL::Bbox_2 simplified_bbox)
{
  // Return true if and only if the non-simplified polygon contains all
  // vertices in the simplified polygon and the order of the vertices in
  // both polygons match. We pass the bounding box of the simplified bounding
  // box as a parameter so that we do not need to construct it in the function
  // body.
  if (non_simplified_polygon.size() < simplified_polygon.size()) {
    return false;  // Simplification cannot create additional vertices
  }

  // If the two bounding boxes neither touch nor overlap, then the
  // non-simplified polygon does not contain the vertices of the simplified
  // polygon. For the if-condition, see:
  // https://stackoverflow.com/questions/325933/determine-whether-two-date-
  // ranges-overlap/325964#325964 (accessed on 2021-10-05).
  CGAL::Bbox_2 non_simplified_bbox = non_simplified_polygon.bbox();
  if (non_simplified_bbox.xmax() < simplified_bbox.xmin()
      || non_simplified_bbox.xmin() > simplified_bbox.xmax()
      || non_simplified_bbox.ymax() < simplified_bbox.ymin()
      || non_simplified_bbox.ymin() > simplified_bbox.ymax()) {
    return false;
  }
  std::vector<unsigned int> indices;
  for (const auto simplified_pt : simplified_polygon) {
    const auto non_simplified_iterator =
      std::find(non_simplified_polygon.vertices_begin(),
                non_simplified_polygon.vertices_end(),
                simplified_pt);

    // Return false if there is no matching vertex in polygon1
    if (non_simplified_iterator == non_simplified_polygon.vertices_end()) {
      return false;
    }
    indices.push_back(distance(non_simplified_polygon.vertices_begin(),
                               non_simplified_iterator));
    if (!std::is_sorted(indices.begin(), indices.end())) {
      return false;
    }
  }
  return std::is_sorted(indices.begin(), indices.end());
}

int simplified_polygon_index(Polygon non_simplified_polygon,
                             std::vector<Polygon> *simplified_polygons,
                             std::vector<CGAL::Bbox_2> *simplified_bboxes,
                             std::list<unsigned int> *unmatched)
{
  // Return index of polygon in simplified_polygons that corresponds to a
  // non-simplified polygon. We pass the bounding boxes of the simplified
  // polygons as a parameter so that the bounding boxes do not need to be
  // calculated when calling contains_vertices_in_order(). We maintain a list
  // of hitherto unmatched polygon indices to avoid unnecessary calls of
  // contains_vertices_in_order().
  for (auto const i : *unmatched) {
    if (contains_vertices_in_order(non_simplified_polygon,
                                   simplified_polygons->at(i),
                                   simplified_bboxes->at(i))) {
      unmatched->remove(i);
      return i;
    }
  }
  return no_matching_non_simplified_polygon;
}

void simplify_inset(InsetState *inset_state)
{
  std::cerr << "Simplifying the inset. "
            << inset_state->n_points()
            << " points in the inset before simplification."
            << std::endl;

  // Store Polygons as CT (Constrained Triangulation) objects. Code inspired
  // by https://doc.cgal.org/latest/Polyline_simplification_2/index.html
  CT ct;
  for (const auto gd : inset_state->geo_divs()) {
    for (auto pwh : gd.polygons_with_holes()) {
      const Polygon &ext_ring = pwh.outer_boundary();
      ct.insert_constraint(ext_ring);
      for (auto it = pwh.holes_begin(); it != pwh.holes_end(); ++it) {
        const Polygon &hole = *it;
        ct.insert_constraint(hole);
      }
    }
  }

  // Simplify polygons
  PS::simplify(ct, Cost(), Stop(0.1));

  std::cerr << "Now doing the bookkeeping" << std::endl;

  // Store each constraint in ct as a polygon. Also store bounding box so
  // that we can match non-simplified and simplified polygons more quickly.
  std::vector<Polygon> simplified_polygons;
  std::vector<CGAL::Bbox_2> simplified_bboxes;
  for (auto it = ct.constraints_begin(); it != ct.constraints_end(); ++it) {

    // First and last point in the constraint are identical. We remove the
    // last point to make the polygon simple.
    Polygon ct_as_polygon(ct.points_in_constraint_begin(*it),
                          --ct.points_in_constraint_end(*it));
    simplified_polygons.push_back(ct_as_polygon);
    simplified_bboxes.push_back(ct_as_polygon.bbox());
  }

  // Match non-simplified polygon to its simplified counterpart
  std::list<unsigned int> unmatched(
    boost::counting_iterator<unsigned int>(0U),
    boost::counting_iterator<unsigned int>(simplified_polygons.size()));
  std::vector<int> matching_simplified_polygon;
  for (const auto gd : inset_state->geo_divs()) {
    for (auto pwh : gd.polygons_with_holes()) {
      Polygon &ext_ring = pwh.outer_boundary();
      const int ext_index = simplified_polygon_index(ext_ring,
                                                     &simplified_polygons,
                                                     &simplified_bboxes,
                                                     &unmatched);
      matching_simplified_polygon.push_back(ext_index);
      for (auto it = pwh.holes_begin(); it != pwh.holes_end(); ++it) {
        const Polygon &hole = *it;
        const int hole_index = simplified_polygon_index(hole,
                                                        &simplified_polygons,
                                                        &simplified_bboxes,
                                                        &unmatched);
        matching_simplified_polygon.push_back(hole_index);
      }
    }
  }

  // Sanity check
  if (std::find(matching_simplified_polygon.begin(),
                matching_simplified_polygon.end(),
                no_matching_non_simplified_polygon) !=
      matching_simplified_polygon.end() ||
      !unmatched.empty()) {
    std::cerr << "Unmatched polygon in simplify_map" << std::endl;
    exit(1);
  }

  unsigned long n_pts = 0;
  for (const auto poly : simplified_polygons) {
    n_pts += poly.size();
  }
  std::cerr << n_pts << " points in ct" << std::endl;
  exit(1);
}
