// TODO: What happens if two polygons have touching lines but the corner
// points are not identical in both lines?

/******************************** Inclusions. ********************************/

// #include <iostream>
// #include <fstream>
// #include <string>
// #include <vector>
// #include <deque>
// #include <map>
// #include <unordered_set>
// #include <unordered_map>
 #include <algorithm>
// #include <chrono>
// #include <cmath>
// #include <iomanip>
//
# include <cctype>
//
// #include "cgal_typedef.h"
#include "inset_state.h"
// #include "geo_div.h"
// #include "polyline_advanced.h"
//
// /**************** 1. Repeat the first point as the last point. ***************/
//
// void repeat_first_point_as_last_point(std::vector<GeoDiv> &gd_vector)
// {
//   int num_vertices = 0;
//   int num_holes = 0;
//   for (GeoDiv &gd : gd_vector) {
//     for (Polygon_with_holes &pgnwh : *gd.ref_to_polygons_with_holes()) {
//       Polygon *outer_pgn = &pgnwh.outer_boundary();
//       outer_pgn->push_back((*outer_pgn)[0]);
//       num_vertices += outer_pgn->size();
//       std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
//       for (auto hole = pgnwh.holes_begin();
//            hole != pgnwh.holes_end();
//            hole++) {
//         hole->push_back((*hole)[0]);
//         num_vertices += hole->size();
//         num_holes += 1;
//       }
//     }
//   }
//   std::cout << std::endl;
//   std::cout << "Original number of vertices (double-counts shared ";
//   std::cout << "boundaries): " << num_vertices << std::endl;
//   std::cout << "Original number of holes: " << num_holes << std::endl;
//   std::cout << std::endl;
//   return;
// }
//
// bool bboxes_overlap_or_touch(Polygon pgn1, Polygon pgn2)
// {
//   CGAL::Bbox_2 bbox1 = pgn1.bbox();
//   CGAL::Bbox_2 bbox2 = pgn2.bbox();
//
//   // Formula to check if two sets of bboxes overlap. See
//   // https://stackoverflow.com/questions/325933/determine-whether-two-date-
//   // ranges-overlap/325964#325964 (accessed on 2021-10-05).
//   return bbox1.xmax() >= bbox2.xmin() && bbox2.xmax() >= bbox1.xmin()
//          && bbox1.ymax() >= bbox2.ymin() && bbox2.ymax() >= bbox1.ymin();
// }
//
// bool pgns_are_separated_and_not_identical(Polygon pgn1, Polygon pgn2)
// {
//   for (Point pt : pgn2) {
//
//     // Test for equality: two polygons are equal if and only if there exists
//     // a cyclic permutation of the vertices of pgn2 such that they are equal
//     // to the vertices of pgn1. See https://doc.cgal.org/latest/Polygon/
//     // classCGAL_1_1Polygon__2.html#a0defc827866d4985ea6f4e6335f8c9a0
//     // (accessed on 2021-10-05).
//     bool identical = (pgn1 == pgn2);
//     auto bounded_side = CGAL::bounded_side_2(pgn1.begin(), pgn1.end(), pt);
//     bool inside_of_pgn = bounded_side == CGAL::ON_BOUNDED_SIDE;
//     bool on_boundary_of_pgn = bounded_side == CGAL::ON_BOUNDARY;
//     if (!identical && (inside_of_pgn || on_boundary_of_pgn)) return false;
//   }
//   return true;
// }
//
// bool check_for_neighbours(std::vector<GeoDiv> gd_vector, Polygon pgn)
// {
//   /**
//    * TODO add set to keep track of polygons (outer/hole)
//    * already registered as part of another polygon
//    */
//
//   // TODO: "continue" is hardly ever a good strategy. Check whether we can
//   // clean it up.
//   for (GeoDiv gd : gd_vector) {
//     for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
//       Polygon outer_pgn = pgnwh.outer_boundary();
//
//       /* If bboxes don't overlap, pgns are not neighbours. Continue. */
//       if (!bboxes_overlap_or_touch(pgn, outer_pgn)) continue;
//
//       /* If pgns are contiguous, pgn is not an island. */
//       if (!pgns_are_separated_and_not_identical(pgn, outer_pgn)) return false;
//
//       std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
//       for (Polygon hole : holes_v) {
//
//         /* If bboxes don't overlap, pgns are not neighbours. Continue. */
//         if (!bboxes_overlap_or_touch(pgn, hole)) continue;
//
//         /* If pgns are contiguous, pgn is not an island. */
//         if (!pgns_are_separated_and_not_identical(pgn, hole)) return false;
//       }
//     }
//   }
//
//   /* If pgn's bbox no overlap and not contiguous, then pgn is an island. */
//   return true;
// }
//
// /********* 2. Return a vector<vector<bool>> of pgns IDed as islands. *********/
//
// std::vector<std::vector<bool> > pgn_is_island(std::vector<GeoDiv> gd_vector)
// {
//
//   std::cout << "Identifying islands..." << std::endl;
//
//   std::vector<std::vector<bool> > pgn_bool_island(gd_vector.size(),
//                                                   std::vector<bool>());
//
//   /* Iterate through all gds and pgnwhs. */
//   for (unsigned int gd_num = 0; gd_num < gd_vector.size(); gd_num++) {
//     GeoDiv gd = gd_vector[gd_num];
//
//     std::cout << "gd: " << gd_num << std::endl;
//
//     std::vector<bool> v_bool(gd.polygons_with_holes().size(), false);
//     pgn_bool_island[gd_num] = v_bool;
//
//     for (unsigned int pgnwh_num = 0;
//          pgnwh_num < gd.polygons_with_holes().size();
//          pgnwh_num++) {
//       Polygon_with_holes pgnwh = gd.polygons_with_holes()[pgnwh_num];
//       Polygon outer_pgn = pgnwh.outer_boundary();
//
//       /* Check if this pgnwh has neighbours. */
//       pgn_bool_island[gd_num][pgnwh_num] =
//         check_for_neighbours(gd_vector, outer_pgn);
//     }
//   }
//
//   /* Count and print total number of islands and non-islands. */
//   int num_islands = 0;
//   int num_non_islands = 0;
//   for (int i = 0; i < (int) pgn_bool_island.size(); i++) {
//     for (int j = 0; j < (int) pgn_bool_island[i].size(); j++) {
//       num_islands = pgn_bool_island[i][j] ? num_islands + 1 : num_islands;
//       num_non_islands = !pgn_bool_island[i][j] ? num_non_islands + 1
//                         : num_non_islands;
//     }
//   }
//   std::cout << "Number of islands: " << num_islands << std::endl;
//   std::cout << "Number of non-islands: " << num_non_islands << std::endl;
//   return pgn_bool_island;
// }
//
// void insert_pll_into_graph(Polyline &pll,
//                            Graph &graph,
//                            Point_vertex_map &pv_map)
// {
//   for (unsigned int i = 0; i < pll.size(); ++i) {
//     vertex_descriptor u, v;
//
//     /* If the point is not yet in the graph. */
//
//     // TODO: pv_map.contains(pll[i]) should work and would be more elegant.
//     // See https://en.cppreference.com/w/cpp/container/map/contains
//     if (pv_map.find(pll[i]) == pv_map.end()) {
//       v = add_vertex(graph);
//       pv_map[pll[i]] = v;
//     } else {
//       v = pv_map[pll[i]];
//     }
//
//     /* Associate the point to the vertex. */
//     graph[v] = pll[i];
//
//     // TODO: This if-condition would better be handled outside the for-loop.
//     // Instead, the loop should start at 1.
//     if (i != 0) {
//       add_edge(u, v, graph);
//     }
//     u = v;
//   }
// }
//
// /*********** 3. Create graph and split graph into unique polylines. **********/
//
// Graph create_pll_graph(std::vector<GeoDiv> gd_vector)
// {
//   Graph graph;
//   Point_vertex_map pv_map;
//   for (GeoDiv gd : gd_vector) {
//     for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
//       Polyline pll_outer;
//       Polygon outer_pgn = pgnwh.outer_boundary();
//       for (Point pt_outer : outer_pgn) {
//         pll_outer.push_back(pt_outer);
//       }
//
//       insert_pll_into_graph(pll_outer, graph, pv_map);
//
//       std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
//       for (Polygon hole : holes_v) {
//         Polyline pll_hole;
//         for (Point pt_hole : hole) {
//           pll_hole.push_back(pt_hole);
//         }
//
//         insert_pll_into_graph(pll_hole, graph, pv_map);
//       }
//     }
//   }
//   return graph;
// }
//
// /* Struct template using Graph internally. */
// template <typename Graph>
//
// // TODO: WE USE UPPER CAMEL CASE FOR STRUCTS: PolylineVisitor
// struct Polyline_visitor
// {
//   std::list<Polyline>& pll_list;
//   const Graph& points_pmap;
//
//   /* Member initializer list. */
//   Polyline_visitor(std::list<Polyline>& lines,
//                    const Graph& points_property_map)
//     : pll_list(lines), points_pmap(points_property_map)
//   {
//   }
//
//   void start_new_polyline() {
//     Polyline V;
//     pll_list.push_back(V);
//   }
//
//   void add_node(typename boost::graph_traits<Graph>::vertex_descriptor vd) {
//
//     // TODO: WOULD IT BE POSSIBLE TO SIMPLIFY THE NEXT TWO LINES TO:
//     // pll_list.push_back(points_pmap[vd]);
//     Polyline& pll = pll_list.back();
//     pll.push_back(points_pmap[vd]);
//   }
//
//   // TODO: COULD WE DO SOMETHING MEANINGFUL IN end_polyline()? FOR EXAMPLE,
//   // DETERMINE WHETHER THE END POINT IS EQUAL TO THE STARTING POINT AND, HENCE
//   // THE POLYLINE COMES FROM AN ISLAND OR A HOLE?
//   void end_polyline() {
//   }
// };
//
// void print_pll_adv(Polyline_advanced pll_adv) {
//   std::cout << pll_adv.gd();
//   std::cout << " | " << pll_adv.pgnwh();
//   std::cout << " | " << pll_adv.pos();
//   std::cout << " | " << pll_adv.is_hole();
//   std::cout << " | " << pll_adv.v1();
//   std::cout << " | " << pll_adv.vl();
//   std::cout << " | " << pll_adv.v2() << std::endl;
// }
//
// /**
//  * Track the matching progress of its polylines to its original gds and pgnwhs.
//  */
// class PgnProg {
// private:
//
// // Can endpts be a simple fixed-size array "Point enpts[2]"?
// std::vector<Point> endpts;
// int island = false;
// int num_holes = 0;
// bool endpts_prog = false;
// bool holes_prog = false;
//
// public:
// void add_endpts(Point endpt1, Point endpt2) {
//   endpts.push_back(endpt1);
//   endpts.push_back(endpt2);
// }
//
// std::vector<Point> get_endpts() {
//   return endpts;
// }
//
// void add_holes(int num_holes_) {
//   num_holes += num_holes_;
// }
//
// void rem_hole() {
//   num_holes--;
// }
//
// void update_prog() {
//   if (num_holes == 0) {
//     holes_prog = true;
//   }
//
//   std::sort(endpts.begin(), endpts.end());
//   int pairs = 0;
//   for (int i = 0; i < (int) endpts.size() - 1; i += 2) {
//     if (endpts[i] == endpts[i + 1]) {
//       pairs++;
//     }
//   }
//
//   /**
//    * If the number of equal endpt pairs is half of the number of endpts,
//    * then the pgnwh has found all its endpt and completed its progress.
//    */
//   if (pairs * 2 == (int) endpts.size()) {
//     endpts_prog = true;
//   }
// }
//
// /**
//    // TODO: BRUSSELS IS NOT ENCLOSED BY A HOLE?
//  * lone_pgn includes islands and polygons enclosed by a hole.
//  * e.g. Brussels
//  *
//    // TODO: I DO NOT UNDERSTAND THE NEXT COMMENT
//  * So, if the current pgnwh has 0 holes, then the pll that it is
//  * matched to is an island, otherwise, it is not an island.
//  */
//
// // TODO: WHY DO WE NEED THIS FUNCTION? IT HARDLY DOES ANYTHING.
// void set_lone_pgn() {
//   island = (num_holes == 0);
// }
//
// /**
//  * To complete its progress, the pgnwh must be an island or it must
//  * fulfill both endpts_prog and holes_prog.
//  */
// bool check_prog() {
//   return island || (endpts_prog && holes_prog);
// }
// };
//
// /**
//  * Create a map where each pgnwh has an associated PgnProg object to track
//  * the matching progress of its polylines to its original gds and pgnwhs.
//  */
// std::map<int, std::map<int, PgnProg> > create_pgnprog_map(
//   std::vector<GeoDiv> gd_vector)
// {
//   std::map<int, std::map<int, PgnProg> > pgnprog_map;
//   for (int i = 0; i < (int) gd_vector.size(); i++) {
//     for (int j = 0; j < (int) gd_vector[i].polygons_with_holes().size(); j++) {
//       Polygon_with_holes pgnwh = gd_vector[i].polygons_with_holes()[j];
//       std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
//       pgnprog_map[i][j].add_holes(holes_v.size());
//     }
//   }
//   return pgnprog_map;
// }
//
// /**
//  * Check if pll_adv is on the boundary of pgn and return the
//  * number of pll_adv vertices on the boundary of the pgn.
//  */
// // TODO: THE FUNCTION NAME DOES NOT REVEAL WHAT THE RETURN VALUE IS. THE
// // FUNCTION RETURNS WHETHER A POLYLINE HAS 0, 1, 2 OR >=3 POINTS ON THE
// // BOUNDARY OF A GIVEN POLYGON.
// int check_if_pll_adv_on_pgn_boundary(Polyline_advanced pll_adv, Polygon pgn)
// {
//   /* Iterate through all vertices (points) inside the polyline. */
//   int num_v_on_boundary = 0;
//   for (int i = 0; i < (int) pll_adv.pll().size(); i++) {
//     Point pt = pll_adv.pll()[i];
//     bool v_on_boundary = CGAL::bounded_side_2(
//       pgn.begin(), pgn.end(), pt) == CGAL::ON_BOUNDARY;
//
//     /* If the vertex is on the pgn's boundary, increment. */
//     if (v_on_boundary) num_v_on_boundary++;
//
//     /**
//      * If the number of iterations has reached 3, break. This means that the
//      * number of vertices on the pgn's boundary may or may not have reached 3.
//      */
//     if (i == 2) break;
//   }
//
//   return num_v_on_boundary;
// }
//
// /**
//  * Check the matching progress of its polylines to its original gds and pgnwhs.
//  */
// void check_pgnprog_map(Polyline_advanced pll_adv,
//                        std::map<int, std::map<int, PgnProg> > &pgnprog_map,
//                        std::vector<int> &plls_match,
//                        bool pgnwh_is_island)
// {
//   int pos = pll_adv.pos();
//   bool is_hole = pll_adv.is_hole();
//
//   /**
//    * If the pll_adv's 1st vertex equals its last vertex and it's a hole,
//    * remove 1 hole from pgnprog_map and increment plls_match by 1.
//    */
//   if (pll_adv.v1() == pll_adv.vl() && is_hole) {
//     pgnprog_map[pll_adv.gd()][pll_adv.pgnwh()].rem_hole();
//     plls_match[pos] += 1;
//
//     /**
//      * If the pll_adv's 1st vertex equals its last vertex and it's an island,
//      * set it to be a lone polygon and increment plls_match by 2.
//      */
//   } else if (pll_adv.v1() == pll_adv.vl() && pgnwh_is_island) {
//     pgnprog_map[pll_adv.gd()][pll_adv.pgnwh()].set_lone_pgn();
//     plls_match[pos] += 2;
//
//     /**
//      * If the pll_adv's 1st vertex equals its last vertex and it's not a hole,
//      * set it to be a lone polygon and increment plls_match by 1.
//      */
//   } else if (pll_adv.v1() == pll_adv.vl() && !is_hole) {
//     pgnprog_map[pll_adv.gd()][pll_adv.pgnwh()].set_lone_pgn();
//     plls_match[pos] += 1;
//
//     /**
//      * Else, add its endpoints, update its progress, and increment plls_match by 1.
//      */
//   } else {
//     pgnprog_map[pll_adv.gd()][pll_adv.pgnwh()]
//     .add_endpts(pll_adv.v1(), pll_adv.vl());
//     pgnprog_map[pll_adv.gd()][pll_adv.pgnwh()].update_prog();
//     plls_match[pos] += 1;
//   }
// }
//
// /**
//  * Depending on num_v_pll_outer, call check_pgnprog_map and plls_adv_by_pos.
//  */
// void check_prog_store_pll(int num_v_on_boundary,
//                           Polyline_advanced pll_adv,
//                           Polygon pgn,
//                           std::map<int, std::vector<Polyline_advanced> > &plls_adv_by_pos,
//                           std::map<int, std::map<int, PgnProg> > &pgnprog_map,
//                           std::vector<int> &plls_match,
//                           bool pgnwh_is_island)
// {
//   /**
//    * If the number of vertices on the pgn's boundary reached 3:
//    * a) Check pgnprog for this pgn taking into account this pll_adv.
//    * b) Store this pll_adv according to its original inset_state position
//    *   along with its associated GeoDiv and Polygon_with_hole.
//    */
//   if (num_v_on_boundary >= 3) {
//     check_pgnprog_map(pll_adv, pgnprog_map, plls_match, pgnwh_is_island);
//     plls_adv_by_pos[pll_adv.pos()].push_back(pll_adv);
//
//     /**
//      * If the number of vertices on the pgn's boundary only reached 2, iterate
//      * through all vertices of the pgn and determine if the 2 common vertices
//      * in the pll_adv and pgn are consecutive, checking both the clockwise and
//      * counter-clockwise orientation of the pll_adv and pgn.
//      */
//   } else if (num_v_on_boundary == 2) {
//     for (int i = 0; i < (int) pgn.size() - 1; i++) {
//       bool direction_1 = pgn[i] == pll_adv.v1() && pgn[i + 1] == pll_adv.v2();
//       bool direction_2 = pgn[i + 1] == pll_adv.v1() && pgn[i] == pll_adv.v2();
//
//       /** If the 2 common vertices are consecutive:
//        * a) Check pgnprog for this pgn taking into account this pll_adv.
//        * b) Store this pll_adv according to its original inset_state position
//        *    along with its associated GeoDiv and Polygon_with_hole.
//        * c) Break because we no longer need to check all other vertices.
//        */
//       if (direction_1 || direction_2) {
//         check_pgnprog_map(pll_adv, pgnprog_map, plls_match, pgnwh_is_island);
//         plls_adv_by_pos[pll_adv.pos()].push_back(pll_adv);
//         break;
//       }
//     }
//   }
// }
//
// /** 6. Match and store polylines according to their original inset_state      **/
// /**    positions along with their associated GeoDiv and Polygon_with_hole.  **/
//
// std::map<int, std::vector<Polyline_advanced> > store_by_pos(
//   std::vector<Polyline> &ct_polylines,
//   std::vector<GeoDiv> gd_vector,
//   std::vector<std::vector<bool> > gd_pgnwh_island_bool)
// {
//   /**
//    * Create map to match and store polylines according to their original
//    * inset_state positions along with their associated GeoDiv and
//    * Polygon_with_hole.
//    */
//   std::map<int, std::vector<Polyline_advanced> > plls_adv_by_pos;
//
//   /**
//    * Create a map where each pgnwh has an associated PgnProg object to track
//    * the matching progress of its polylines to its original gds and pgnwhs.
//    */
//   std::map<int, std::map<int, PgnProg> >
//   pgnprog_map = create_pgnprog_map(gd_vector);
//
//   /**
//    * Create a vector of polylines and an associated value x,
//    * where x=1 or x=2, and x represents the number of polygons
//    * it has been matched to.
//    */
//   std::vector<int> plls_match(ct_polylines.size(), 0);
//
//   std::cout << "Matching polylines to polygons..." << std::endl;
//
//   /* Iterate through all gds and pgnwhs and then through all polylines. */
//   for (int gd_num = 0; gd_num < (int) gd_vector.size(); gd_num++) {
//     std::cout << "gd: " << gd_num << std::endl;
//     for (int pgnwh_num = 0;
//          pgnwh_num < (int) gd_vector[gd_num].polygons_with_holes().size();
//          pgnwh_num++) {
//
//       for (int pos = 0; pos < (int) ct_polylines.size(); pos++) {
//         Polyline pll = ct_polylines[pos];
//
//         /* If the pgnwh has already found/matched all its polylines, break */
//         if (pgnprog_map[gd_num][pgnwh_num].check_prog()) break;
//
//         /* If the pll has been fully matched, continue */
//         if (plls_match[pos] == 2) continue;
//
//         Polygon_with_holes pgnwh = gd_vector[gd_num]
//                                    .polygons_with_holes()[pgnwh_num];
//
//         bool pgnwh_is_island = gd_pgnwh_island_bool[gd_num][pgnwh_num];
//
//         /**
//          * Create a Polyline_advanced object with the current polyline's position,
//          * Polyline object, gd_num, pgnwh_num, and whether the pll
//          * is a hole or not.
//          */
//         Polyline_advanced pll_adv_outer(pos, pll, gd_num, pgnwh_num, false);
//         Polygon outer_pgn = pgnwh.outer_boundary();
//
//         /**
//          * Check if pll_adv_outer is on the boundary of outer_pgn and return
//          * the number of pll vertices on the boundary of the pgn.
//          */
//         int num_v_pll_adv_outer = check_if_pll_adv_on_pgn_boundary(pll_adv_outer,
//                                                                    outer_pgn);
//
//         /**
//          * Depending on num_v_pll_adv_outer, call check_pgnprog_map and plls_adv_by_pos.
//          */
//         check_prog_store_pll(num_v_pll_adv_outer,
//                              pll_adv_outer,
//                              outer_pgn,
//                              plls_adv_by_pos,
//                              pgnprog_map,
//                              plls_match,
//                              pgnwh_is_island);
//
//         std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
//         for (Polygon hole : holes_v) {
//           Polyline_advanced pll_adv_hole(pos, pll, gd_num, pgnwh_num, true);
//
//           /**
//            * Check if pll_adv_hole is on the boundary of hole and return
//            * the number of pll vertices on the boundary of the pgn.
//            */
//           int num_v_pll_adv_hole = check_if_pll_adv_on_pgn_boundary(pll_adv_hole, hole);
//
//           /**
//            * Depending on num_v_pll_adv_hole, call check_pgnprog_map and plls_adv_by_pos.
//            */
//           check_prog_store_pll(num_v_pll_adv_hole,
//                                pll_adv_hole,
//                                hole,
//                                plls_adv_by_pos,
//                                pgnprog_map,
//                                plls_match,
//                                pgnwh_is_island);
//         }
//       }
//     }
//   }
//   std::cout << std::endl;
//
//   return plls_adv_by_pos;
// }
//
// /** 8. Store polylines according to their GeoDivs and Polygon_with_holes    **/
// /**    along with their associated ct_polyline positions.                   **/
//
// std::map<int, std::map<int, std::vector<Polyline_advanced> > >
// store_by_gd_pgnwh(std::vector<GeoDiv> gd_vector,
//                   CT &ct, std::map<int,
//                                    std::vector<Polyline_advanced> > &plls_adv_by_pos)
// {
//   std::map<int, std::map<int, std::vector<Polyline_advanced> > > plls_adv_by_gd_pgnwh;
//
//   for (int gd_num = 0; gd_num < (int) gd_vector.size(); gd_num++) {
//     for (int pgnwh_num = 0;
//          pgnwh_num < (int) gd_vector[gd_num].polygons_with_holes().size();
//          pgnwh_num++) {
//       /**
//        * Iterates through each cit, which represents a now-simplified polyline.
//        */
//       int cit_num = 0;
//       for (auto cit = ct.constraints_begin(); cit != ct.constraints_end();
//            cit++) {
//         /* Iterates through each non-simplified polyline. */
//         for (Polyline_advanced pll_adv : plls_adv_by_pos[cit_num]) {
//
//           /* If this polyline's associated gd and pgnwh are valid: */
//           if (pll_adv.gd() == gd_num && pll_adv.pgnwh() == pgnwh_num) {
//             Polyline pll_ct;
//             for (auto v : ct.vertices_in_constraint(*cit)) {
//               pll_ct.push_back(v->point());
//             }
//             Polyline_advanced pll_adv_new(pll_adv.pos(),
//                                           pll_ct,
//                                           pll_adv.gd(),
//                                           pll_adv.pgnwh(),
//                                           pll_adv.is_hole());
//
//             /* Stores polyline in plls_adv_by_gd_pgnwh. */
//             plls_adv_by_gd_pgnwh[gd_num][pgnwh_num].push_back(pll_adv_new);
//             break;
//           }
//         }
//         cit_num++;
//       }
//     }
//   }
//   return plls_adv_by_gd_pgnwh;
// }
//
// /******* 9. Set all polylines to not-visited and sort pll_adv_by_gd_pgnwh. *******/
//
// void set_visited_vals(std::unordered_map<int,
//                                          std::unordered_map<int, std::unordered_map<int, bool> > > &visited,
//                       std::map<int, std::map<int, std::vector<Polyline_advanced> > > &plls_adv_by_gd_pgnwh)
// {
//   for (auto [gd_num, map_iv] : plls_adv_by_gd_pgnwh) {
//     for (auto [pgnwh_num, pll_adv_v] : map_iv) {
//       for (Polyline_advanced pll_adv : pll_adv_v) {
//         /* Set all polylines to not-visited. */
//         visited[gd_num][pgnwh_num][pll_adv.pos()] = false;
//       }
//
//       /**
//        * Sort each pll_adv_v so that all holes are in front
//        * to evaluate them in the next step first.
//        */
//       std::sort(plls_adv_by_gd_pgnwh[gd_num][pgnwh_num].begin(),
//                 plls_adv_by_gd_pgnwh[gd_num][pgnwh_num].end(),
//                 [](Polyline_advanced pll_adv_1, Polyline_advanced pll_adv_2) {
//         return pll_adv_1.is_hole() > pll_adv_2.is_hole();
//       });
//     }
//   }
// }
//
// /**
//  * Check if the holes are not inside the pgn by checking if
//  * holes' middle vertices are not inside the pgn's boundary.
//  *
//  * Then update gd_final accordingly.
//  */
//
// void check_holes_inside_pgn(std::vector<Polygon> &holes_v,
//                             Polygon pgn,
//                             GeoDiv &gd_final)
// {
//   bool holes_inside = true;
//   for (Polygon hole : holes_v) {
//     if (CGAL::bounded_side_2(pgn.begin(), pgn.end(), hole[hole.size() / 2])
//         != CGAL::ON_BOUNDED_SIDE) {
//       holes_inside = false;
//     }
//   }
//
//   if (holes_inside) {
//     Polygon_with_holes pgnwh_final(pgn, holes_v.begin(), holes_v.end());
//
//     /* Update gd_final by adding pgnwh. */
//     gd_final.push_back(pgnwh_final);
//
//     /* No other pgnwh would have the same holes. */
//     holes_v.clear();
//   }
// }
//
// /* Connect together all polylines belonging to a pgnwh. */
//
// void connect_polylines_in_deq(
//   std::map<int, std::map<int, std::vector<Polyline_advanced> > > plls_adv_by_gd_pgnwh,
//   std::unordered_map
//   <int, std::unordered_map<int, std::unordered_map<int, bool> > > &visited,
//   std::deque<Polyline_advanced> &pll_adv_deq,
//   Polyline_advanced pll_adv,
//   int gd_num,
//   int pgnwh_num)
// {
//   /* Search for polylines that connect with the plls inside the pll_adv_deq. */
//   while (1) {
//     bool pll_adv_siblings_found = false;
//
//     /* Iterate through all plls inside the current gd/pgnwh. */
//     for (Polyline_advanced pll_adv_inner : plls_adv_by_gd_pgnwh[pll_adv.gd()][pll_adv.pgnwh()]) {
//       if (visited[gd_num][pgnwh_num][pll_adv_inner.pos()] == true) {
//         continue;
//       }
//
//       /**
//        * If pll_adv_deq's 1st Polyline_advanced's 1st vertex == pll_adv_inner's last vertex,
//        * connect them.
//        */
//       if (pll_adv_deq.front().v1()[0] == pll_adv_inner.vl()[0] &&
//           pll_adv_deq.front().v1()[1] == pll_adv_inner.vl()[1]) {
//         pll_adv_deq.front().pop_front(); // Remove otherwise repeated vertex
//         pll_adv_deq.push_front(pll_adv_inner);
//         visited[gd_num][pgnwh_num][pll_adv_inner.pos()] = true;
//         pll_adv_siblings_found = true;
//
//         /**
//          * If pll_adv_deq's last Polyline_advanced's last vertex == pll_adv_inner's 1st vertex,
//          * connect them.
//          */
//       } else if (pll_adv_deq.back().vl()[0] == pll_adv_inner.v1()[0] &&
//                  pll_adv_deq.back().vl()[1] == pll_adv_inner.v1()[1]) {
//         pll_adv_deq.back().pop_back(); // Remove otherwise repeated vertex
//         pll_adv_deq.push_back(pll_adv_inner);
//         visited[gd_num][pgnwh_num][pll_adv_inner.pos()] = true;
//         pll_adv_siblings_found = true;
//
//         /**
//          * If pll_adv_deq's 1st Polyline_advanced's 1st vertex == pll_adv_inner's 1st vertex,
//          * reverse pll_adv_inner and connect them.
//          */
//       } else if (pll_adv_deq.front().v1()[0] == pll_adv_inner.v1()[0] &&
//                  pll_adv_deq.front().v1()[1] == pll_adv_inner.v1()[1]) {
//         Polyline pll_new = pll_adv_inner.pll();
//         std::reverse(pll_new.begin(), pll_new.end());
//         pll_adv_inner.set_pll(pll_new);
//         pll_adv_deq.front().pop_front(); // Remove otherwise repeated vertex
//         pll_adv_deq.push_front(pll_adv_inner);
//         visited[gd_num][pgnwh_num][pll_adv_inner.pos()] = true;
//         pll_adv_siblings_found = true;
//
//         /**
//          * If pll_adv_deq's last Polyline_advanced's last vertex == pll_adv_inner's last vertex,
//          * reverse pll_adv_inner and connect them.
//          */
//       } else if (pll_adv_deq.back().vl()[0] == pll_adv_inner.vl()[0] &&
//                  pll_adv_deq.back().vl()[1] == pll_adv_inner.vl()[1]) {
//         Polyline pll_new = pll_adv_inner.pll();
//         std::reverse(pll_new.begin(), pll_new.end());
//         pll_adv_inner.set_pll(pll_new);
//         pll_adv_deq.back().pop_back(); // Remove otherwise repeated vertex
//         pll_adv_deq.push_back(pll_adv_inner);
//         visited[gd_num][pgnwh_num][pll_adv_inner.pos()] = true;
//         pll_adv_siblings_found = true;
//       }
//     }
//     /**
//      * If siblings cannot be found for any while iteration, then the current
//      * pll_adv_deq along with its plls have no more polylines to connect to. So, break.
//      */
//     if (!pll_adv_siblings_found) break;
//   }
// }
//
// /******************* 10. Assemble polylines into polygons. *******************/
//
// void assemble_plls_adv_to_pgn(
//   std::map<int, std::map<int, std::vector<Polyline_advanced> > > &plls_adv_by_gd_pgnwh,
//   std::unordered_map
//   <int, std::unordered_map<int, std::unordered_map<int, bool> > > &visited,
//   std::vector<GeoDiv> &gd_vector_final,
//   std::vector<GeoDiv> gd_vector_org)
// {
//   std::cout << "Assembling polylines into polygons..." << std::endl;
//   for (auto [gd_num, map_iv] : plls_adv_by_gd_pgnwh) {
//     GeoDiv gd_final(gd_vector_org[gd_num].id());
//
//     for (auto [pgnwh_num, pll_adv_v] : map_iv) {
//
//       /* 1 holes_v per pgnwh. */
//       std::vector<Polygon> holes_v;
//
//       for (Polyline_advanced pll_adv : pll_adv_v) {
//         Polygon hole_or_island; // This will only be for islands/holes anyway
//
//         if (visited[gd_num][pgnwh_num][pll_adv.pos()]) continue;
//
//         /* If it is a single polyline (e.g. island, hole): */
//         if (pll_adv.v1() == pll_adv.vl()
//             && !visited[gd_num][pgnwh_num][pll_adv.pos()]) {
//
//           visited[gd_num][pgnwh_num][pll_adv.pos()] = true;
//
//           for (Point pt : pll_adv.pll()) {
//             hole_or_island.push_back(pt);
//           }
//
//           /* If the pll_adv is a hole, add it to holes_v. */
//           if (pll_adv.is_hole()) {
//             holes_v.push_back(hole_or_island);
//
//             // TODO
//             // - add special case where the hole contains a hole.
//
//             /* If the pll_adv is not a hole (and so it's an island): */
//           } else {
//
//             /* If there does not exist a hole inside the pgnwh: */
//             if (holes_v.empty()) {
//               Polygon_with_holes pgnwh(hole_or_island);
//               gd_final.push_back(pgnwh);
//
//               /* If there exists a hole inside the pgnwh: */
//             } else {
//               check_holes_inside_pgn(holes_v, hole_or_island, gd_final);
//             }
//           }
//
//           /* If it is part of a pgnwh with >1 polyline: */
//         } else {
//           /* pll_adv_deq represents the outer_boundary of a pgnwh. */
//           std::deque<Polyline_advanced> pll_adv_deq;
//           pll_adv_deq.push_back(pll_adv);
//           visited[gd_num][pgnwh_num][pll_adv.pos()] = true;
//
//           connect_polylines_in_deq(plls_adv_by_gd_pgnwh,
//                                    visited, pll_adv_deq, pll_adv,
//                                    gd_num,
//                                    pgnwh_num);
//
//           /* With the finished deque, create the pgnwh's outer_boundary(). */
//           Polygon outer_pgn;
//           for (Polyline_advanced pll_adv : pll_adv_deq) {
//             for (Point pt : pll_adv.pll()) {
//               outer_pgn.push_back(pt);
//             }
//           }
//
//           /* If the pll_adv is a hole, add it to holes_v. */
//           /* Special case of a hole having >1 polylines e.g. Switzerland */
//           if (pll_adv.is_hole()) {
//             holes_v.push_back(outer_pgn);
//             continue;
//           }
//
//           /* If there does not exist a hole inside the pgnwh: */
//           if (holes_v.empty()) {
//             Polygon_with_holes pgnwh_final(outer_pgn);
//             gd_final.push_back(pgnwh_final);
//
//             /* If there exists a hole inside the pgnwh: */
//           } else {
//             check_holes_inside_pgn(holes_v, outer_pgn, gd_final);
//           }
//         }
//       }
//     }
//     gd_vector_final.push_back(gd_final);
//   }
// }
//
// /* Print the number of points in the gd_vector. */
//
// void print_num_pts(std::vector<GeoDiv> gd_vector) {
//   int num_pts = 0;
//   for (GeoDiv gd : gd_vector) {
//     for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
//       Polygon outer_pgn = pgnwh.outer_boundary();
//       num_pts += outer_pgn.size();
//
//       std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
//       for (Polygon hole : holes_v)
//         num_pts += hole.size();
//     }
//   }
//   std::cout << num_pts << std::endl;
// }
//
// /*********************** 11. Remove the the last point. **********************/
//
// void remove_first_point_as_last_point(std::vector<GeoDiv> &gd_vector) {
//   for (GeoDiv &gd : gd_vector) {
//     for (Polygon_with_holes &pgnwh : *gd.ref_to_polygons_with_holes()) {
//
//       Polygon *outer_pgn = &pgnwh.outer_boundary();
//       auto outer_pgn_it = outer_pgn->end();
//       outer_pgn_it--;
//       outer_pgn->erase(outer_pgn_it);
//
//       for (auto hole = pgnwh.holes_begin(); hole != pgnwh.holes_end(); hole++)
//       {
//         auto hole_it = hole->end();
//         hole_it--;
//         hole->erase(hole_it);
//       }
//     }
//   }
// }
//
// /***************** 12. Orientate exterior and interior rings. ****************/
//
// void orientate_exterior_and_interior_rings(std::vector<GeoDiv> &gd_vector) {
//   for (GeoDiv &gd : gd_vector) {
//     for (Polygon_with_holes &pgnwh : *gd.ref_to_polygons_with_holes()) {
//
//       Polygon *ext_ring = &pgnwh.outer_boundary();
//
//       if (ext_ring->is_clockwise_oriented()) {
//         ext_ring->reverse_orientation();
//       }
//
//       for (auto hole = pgnwh.holes_begin(); hole != pgnwh.holes_end(); hole++) {
//         if (hole->is_counterclockwise_oriented()) {
//           hole->reverse_orientation();
//         }
//       }
//     }
//   }
// }

void simplify_map(InsetState *inset_state)
{
  // TODO: Should steps 1 and 2 be swapped? Currently, the polygons are no
  // longer simple (at least in CGAL terms) when step 2 is entered.
  // Because the function pgns_are_not_contiguous(), called from
  // pgn_is_island(), makes checks about ON_BOUNDED_SIDE, I suspect we can
  // run into undefined behavior.

  /********************************* Steps: **********************************/
  /* 1. Repeat the first point as the last point.                            */
  /* 2. Return a vector<vector<bool>> of pgns IDed as islands.               */
  /* 3. Create graph and split graph into unique polylines.                  */
  /* 4. Store polylines in a CT object from polyline_list.                   */
  /* 5. Store polylines in a vector in the order of CT polylines.            */
  /* 6. Match and store polylines according to their original inset_state    */
  /*    positions along with their associated GeoDiv and Polygon_with_hole.  */
  /* 7. Simplify polylines.                                                  */
  /* 8. Store polylines according to their GeoDivs and Polygon_with_holes    */
  /*    along with their associated original inset_state positions.          */
  /* 9. Set all polylines to not-visited and sort pll_adv_by_gd_pgnwh.       */
  /* 10. Assemble polylines into polygon.                                    */
  /* 11. Remove the last point.                                              */
  /* 12. Orientate all exterior counterclockwise and interior rings          */
  /*     clockwise.                                                          */

  // std::vector<GeoDiv> gd_vector = inset_state->geo_divs();
  //
  // /* 1. Repeat the first point as the last point.                            */
  // repeat_first_point_as_last_point(gd_vector);
  //
  // /* Start timer for step 2. */
  // const auto start_s2 = std::chrono::system_clock::now();
  //
  // /* 2. Return a vector<vector<bool>> of pgns IDed as islands                */
  // std::vector<std::vector<bool> > pgn_bool_island =
  //   pgn_is_island(gd_vector);
  //
  // /* Get and print step 2's elapsed time */
  // const std::chrono::duration<double, std::milli> dur_s2 =
  //   std::chrono::system_clock::now() - start_s2;
  // std::cout << "get_gd_pgnwh_island_bool() time elapsed: ";
  // std::cout << dur_s2.count() << " ms (";
  // std::cout << dur_s2.count() / 1000 << " s)" << std::endl;
  // std::cout << std::endl;
  //
  // /* Start timer for remaining steps. */
  // const auto start_s311 = std::chrono::system_clock::now();
  //
  // /* 3. Create graph and split graph into unique polylines.                  */
  // Graph graph = create_pll_graph(gd_vector);
  // std::list<Polyline> polyline_list;
  // Polyline_visitor<Graph> polyline_visitor(polyline_list, graph);
  // CGAL::split_graph_into_polylines(graph, polyline_visitor);
  //
  // for (auto const &pll : polyline_list) {
  //   std::cerr << "A new polyline:" << std::endl;
  //   for (auto pt : pll) {
  //     std::cerr << "\t" << pt << std::endl;
  //   }
  // }
  //
  // exit(1);

  /* 4. Store polylines from polyline_list in a CT object.                   */
  CT ct;
  // for (Polyline polyline : polyline_list) {
  //   ct.insert_constraint(polyline.begin(), polyline.end());
  // }

  // MTG: Store Polygons as CT objects. Code inspired by
  // https://doc.cgal.org/latest/Polyline_simplification_2/index.html
  for (const auto gd : inset_state->geo_divs()) {
    std::cerr << gd.id() << std::endl;
    for (auto pwh : gd.polygons_with_holes()) {
      //std::cerr << pwh << std::endl;
      Polygon &ext_ring = pwh.outer_boundary();
      ct.insert_constraint(ext_ring);
      for (Polygon_with_holes::Hole_const_iterator it = pwh.holes_begin();
           it != pwh.holes_end();
           ++it) {
        const Polygon &hole = *it;
        ct.insert_constraint(hole);
      }
    }
  }

  /* 5. Store polylines in a vector in the order of CT polylines.            */
  /**
   * This is so that we no longer need to use the CT data structure when
   * carrying out operations, which has limited methods available to use.
   * Instead, we can use the vector data structure and its methods.
   */
//   std::vector<Polyline> ct_polylines;
//   for (auto cit = ct.constraints_begin(); cit != ct.constraints_end()
//        ; cit++) {
//     Polyline pll;
//     for (auto v : ct.vertices_in_constraint(*cit)) {
//       pll.push_back(v->point());
//     }
//     ct_polylines.push_back(pll);
//   }
//
//   /* 6. Match and store polylines according to their original inset_state      */
//   /*    positions along with their associated GeoDiv and Polygon_with_hole.  */
//   std::map<int, std::vector<Polyline_advanced> >
//   plls_adv_by_pos = store_by_pos(ct_polylines, gd_vector, pgn_bool_island);
//
/* 7. Simplify polylines.                                                  */

  // Simplify polygons
  PS::simplify(ct, Cost(), Stop(0.1));

  // Store each constraint in ct as a polygon
  std::vector<Polygon> simplified_polygons;
  for (Constraint_iterator cit = ct.constraints_begin();
       cit != ct.constraints_end();
       ++cit) {

    // First and last point in the constraint are identical. We remove the
    // last point to make the polygon simple.
    Polygon ct_as_polygon(ct.points_in_constraint_begin(*cit),
                          --ct.points_in_constraint_end(*cit));
    simplified_polygons.push_back(ct_as_polygon);
  }
  for (const auto gd : inset_state->geo_divs()) {
    std::cerr << "GeoDiv: " << gd.id() << std::endl;
    for (auto pwh : gd.polygons_with_holes()) {
      Polygon &ext_ring = pwh.outer_boundary();

      std::cerr << "An exterior ring in GeoDiv" << gd.id() << ":" << std::endl;
      for (Polygon::Vertex_iterator vit = ext_ring.vertices_begin();
           vit != ext_ring.vertices_end();
           ++vit) {
        std::cerr << *vit
                  << " ("
                  << (*vit).x() - (int)(*vit).x()
                  << " "
                  << (*vit).y() - (int)(*vit).y()
                  << ")"
                  << std::endl;
      }

      for (std::vector<Polygon>::iterator pit = simplified_polygons.begin();
           pit != simplified_polygons.end();
           ++pit) {
        Polygon simplified_polygon = *pit;
        std::cerr << "A simplified polygon: " << std::endl;
        for (Polygon::Vertex_iterator spit = simplified_polygon.vertices_begin();
             spit != simplified_polygon.vertices_end();
             ++spit) {
          std::cerr << *spit
                    << " ("
                    << (*spit).x() - (int)(*spit).x()
                    << " "
                    << (*spit).y() - (int)(*spit).y()
                    << ")"
                    << std::endl;
        }
        bool simplified_polygon_in_ext_ring =
          std::includes(ext_ring.vertices_begin(),
                        ext_ring.vertices_end(),
                        simplified_polygon.vertices_begin(),
                        simplified_polygon.vertices_end());
        if (simplified_polygon_in_ext_ring) {
          std::cerr << "An exterior ring of "
                    << gd.id()
                    << " includes a CT"
                    << std::endl;
        }
        //       std::cerr << "New constraint\n" << std::endl;
        //       for (Points_in_constraint_iterator vit =
        //              ct.points_in_constraint_begin(*cit);
        //            vit != ct.points_in_constraint_end(*cit);
        //            ++vit) {
        //         std::cerr << *vit << std::endl;
        //       }
      }
    }
  }

  const std::vector<char> v1 = {'x', 'h', 'f', 'b', 'c', 'a'},
                          v2 = {'a', 'b', 'c'},
                          v3 = {'a', 'c'},
                          v4 = {'a', 'a', 'b'},
                          v5 = {'g'},
                          v6 = {'x', 'f', 'a'};
  std::cerr << "v1 includes v6: "
            << std::boolalpha
            << std::includes(v1.begin(), v1.end(), v6.begin(), v6.end())
            << std::endl;

  //std::cerr << v1 << "\nincludes:\n" << std::boolalpha
  //   << v2 << ": " << std::includes(v1.begstd::includes(v1.begin(), v1.end(), v2.begin(), v2.end())in(), v1.end(), v2.begin(), v2.end()) << '\n'
  //   << v3 << ": " << std::includes(v1.begin(), v1.end(), v3.begin(), v3.end()) << '\n'
  //   << v4 << ": " << std::includes(v1.begin(), v1.end(), v4.begin(), v4.end()) << '\n'
  //   << v5 << ": " << std::includes(v1.begin(), v1.end(), v5.begin(), v5.end()) << '\n'
  //   << v6 << ": " << std::includes(v1.begin(), v1.end(), v6.begin(), v6.end()) << '\n'
  //        << std::endl;

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
