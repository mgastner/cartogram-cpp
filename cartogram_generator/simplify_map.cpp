// TODO: What happens if two polygons have touching lines but the corner
// points are not identical in both lines?

/******************************** Inclusions. ********************************/
#include <algorithm>
#include "inset_state.h"



bool contains_vertices_in_order(Polygon polygon1, Polygon polygon2)
{
  // Return true if and only if polygon1 contains all vertices in polygon2 and
  // the order of the vertices in both polygons match
  if (polygon1.size() < polygon2.size()) {
    return false;
  }
  std::vector<unsigned int> indices;
  for (const auto pt2 : polygon2) {
    const auto pt1_iterator = std::find(polygon1.vertices_begin(),
                                        polygon1.vertices_end(),
                                        pt2);

    // Return false if there is no matching vertex in polygon1
    if (pt1_iterator == polygon1.vertices_end()) {
      return false;
    }
    indices.push_back(distance(polygon1.vertices_begin(), pt1_iterator));
    if (!std::is_sorted(indices.begin(), indices.end())) {
      return false;
    }
  }
  return std::is_sorted(indices.begin(), indices.end());
}

int simplified_polygon_index(Polygon non_simplified,
                             std::vector<Polygon> *simplified,
                             std::list<unsigned int> *unmatched)
{
  // Return index of polygon in *simplified that corresponds to a
  // non-simplified polygon
  for (auto const i : *unmatched) {
    if (contains_vertices_in_order(non_simplified, simplified->at(i))) {
      unmatched->remove(i);
      return i;
    }
  }
  return -1;  // -1 implies that no matching simplified polygon was found
}

void simplify_map(InsetState *inset_state)
{

  std::cerr << "In simplfy_map()" << std::endl;

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

  // Store each constraint in ct as a polygon
  std::vector<Polygon> simplified_polygons;
  for (auto it = ct.constraints_begin(); it != ct.constraints_end(); ++it) {

    // First and last point in the constraint are identical. We remove the
    // last point to make the polygon simple.
    Polygon ct_as_polygon(ct.points_in_constraint_begin(*it),
                          --ct.points_in_constraint_end(*it));
    simplified_polygons.push_back(ct_as_polygon);
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
                                                     &unmatched);
      matching_simplified_polygon.push_back(ext_index);
      for (auto it = pwh.holes_begin(); it != pwh.holes_end(); ++it) {
        const Polygon &hole = *it;
        const int hole_index = simplified_polygon_index(hole,
                                                        &simplified_polygons,
                                                        &unmatched);
        matching_simplified_polygon.push_back(hole_index);
      }
    }
  }
  if (std::find(matching_simplified_polygon.begin(),
                matching_simplified_polygon.end(),
                -1) !=
      matching_simplified_polygon.end() ||
      !unmatched.empty()) {
    std::cerr << "Unmatched polygon in simplify_map" << std::endl;
    exit(1);
  }
  exit(1);

//
//   /* 8. Store polylines according to their GeoDivs and Polygon_with_holes    */
//   /*    along with their associated original inset_state positions.          */
//   std::map<int, std::map<int, std::vector<Polyline_advanced> > >
//   plls_adv_by_gd_pgnwh = store_by_gd_pgnwh(gd_vector, ct, plls_adv_by_pos);
//
//   /* 9. Set all polylines to not-visited and sort pll_adv_by_gd_pgnwh.       */
//   std::unordered_map<int,
//                      std::unordered_map<int, std::unordered_map<int, bool> >
//                      > visited;
//   set_visited_vals(visited, plls_adv_by_gd_pgnwh);
//
//   /* 10. Assemble polylines into polygons.                                   */
//   std::vector<GeoDiv> gd_vector_simp;
//   assemble_plls_adv_to_pgn(plls_adv_by_gd_pgnwh, visited, gd_vector_simp, gd_vector);
//
//   /* Print number of points before simplifying (gd_vector). */
//   std::cout << "Number of vertices before simplifying (gd_vector): ";
//   print_num_pts(gd_vector);
//
//   /* Print number of points after simplifying. */
//   std::cout << "Number of vertices after simplifying: ";
//   print_num_pts(gd_vector_simp);
//
//   /* 11. Remove the last point.                                              */
//   remove_first_point_as_last_point(gd_vector_simp);
//
//   /* 12. Orientate all exterior rings counterclockwise and interior rings    */
//   /*     clockwise.                                                          */
//   orientate_exterior_and_interior_rings(gd_vector_simp);
//
//   /* Set gd_vector_simp as inset_state's gd_vector. */
//   inset_state->set_geo_divs(gd_vector_simp);
//
//   const std::chrono::duration<double, std::milli>
//   dur_s311 = std::chrono::system_clock::now() - start_s311;
//   std::cout << "Remaining simplification steps time elapsed: ";
//   std::cout << dur_s311.count() << " ms (";
//   std::cout << dur_s311.count() / 1000 << " s)" << std::endl;
//   std::cout << std::endl;
//
//   double dur_total = dur_s2.count() + dur_s311.count();
//   std::cout << "Total simplification time elapsed: " << dur_total << " ms (";
//   std::cout << dur_total / 1000 << " s)" << std::endl;
//   std::cout << std::endl;
}
