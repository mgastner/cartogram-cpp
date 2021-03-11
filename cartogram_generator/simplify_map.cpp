/******************************** Inclusions. ********************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>

#include "cgal_typedef.h"
#include "map_state.h"
#include "geo_div.h"
#include "pll.h"
#include "densify.h"

/********** 1. Function to repeat the first point as the last point. *********/

void repeat_first_point_as_last_point(std::vector<GeoDiv> &gd_vector)
{
  int num_vertices = 0;
  int num_holes = 0;

  for (GeoDiv &gd : gd_vector) {
    for (Polygon_with_holes &pgnwh : *gd.ref_to_polygons_with_holes()) {

      Polygon *outer_pgn = &pgnwh.outer_boundary();
      outer_pgn->push_back((*outer_pgn)[0]);

      num_vertices += outer_pgn->size();

      std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
      for (auto hole = pgnwh.holes_begin(); hole != pgnwh.holes_end(); hole++) {
        hole->push_back((*hole)[0]); 

        num_vertices += hole->size();
        num_holes += 1;
      }
    }
  }
  std::cout << std::endl;
  std::cout << "Original number of vertices (double-counts shared ";
  std::cout << "boundaries): " << num_vertices << std::endl;
  std::cout << "Original number of holes: " << num_holes << std::endl;
  std::cout << std::endl;
}

bool bboxes_overlap(Polygon pgn1, Polygon pgn2)
{
  CGAL::Bbox_2 bbox1 = pgn1.bbox();
  CGAL::Bbox_2 bbox2 = pgn2.bbox();

  /* Formula to check if two sets of bboxes overlap. */
  return bbox1.xmax() >= bbox2.xmin() && bbox2.xmax() >= bbox1.xmin()
    && bbox1.ymax() >= bbox2.ymin() && bbox2.ymax() >= bbox1.ymin();
}

bool pgns_are_not_contiguous(Polygon pgn1, Polygon pgn2)
{
  for (Point pt : pgn2) {
    bool identical = pgn1 == pgn2;
    auto bounded_side = CGAL::bounded_side_2(pgn1.begin(), pgn1.end(), pt);
    bool inside_of_pgn = bounded_side == CGAL::ON_BOUNDED_SIDE;
    bool on_boundary_of_pgn = bounded_side == CGAL::ON_BOUNDARY;
    if (!identical && (inside_of_pgn || on_boundary_of_pgn)) return false;
  }
  return true;
}

bool check_for_neighbours(std::vector<GeoDiv> gd_vector, Polygon pgn)
{
  /**
   * TODO add set to keep track of polygons (outer/hole) 
   * already registered as part of another polygon
   */

  for (GeoDiv gd : gd_vector) {
    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      Polygon outer_pgn = pgnwh.outer_boundary();

      /* If bboxes don't overlap, pgns are not neighbours. Continue. */
      if (!bboxes_overlap(pgn, outer_pgn)) continue;

      /* If pgns are contiguous, pgn is not an island. */
      if (!pgns_are_not_contiguous(pgn, outer_pgn)) return false;

      std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
      for (Polygon hole : holes_v) {

        /* If bboxes don't overlap, pgns are not neighbours. Continue. */
        if (!bboxes_overlap(pgn, hole)) continue;

        /* If pgns are contiguous, pgn is not an island. */
        if (!pgns_are_not_contiguous(pgn, hole)) return false;
      }
    }
  }

  /* If pgn's bbox no overlap and not contiguous, then pgn is an island. */
  return true;
}

/*** 2. Function to return a vector<vector<bool>> of pgns IDed as islands. ***/

std::vector<std::vector<bool>> get_pgn_bool_island(
    std::vector<GeoDiv> gd_vector)
{

  std::cout << "Identifying islands..." << std::endl;

  std::vector<std::vector<bool>> pgn_bool_island(gd_vector.size(),
                                                 std::vector<bool>());

  /* Iterate through all gds and pgnwhs. */
  for (int gd_num = 0;  gd_num < (int) gd_vector.size(); gd_num++) {
    GeoDiv gd = gd_vector[gd_num];

    std::cout << "gd: " << gd_num << std::endl;

    std::vector<bool> v_bool(gd.polygons_with_holes().size(), false);
    pgn_bool_island[gd_num] = v_bool;

    for (int pgnwh_num = 0; pgnwh_num < (int) gd.polygons_with_holes().size() 
                          ; pgnwh_num++) {

      Polygon_with_holes pgnwh = gd.polygons_with_holes()[pgnwh_num];
      Polygon outer = pgnwh.outer_boundary();

      /* Check if this pgnwh has neighbours. */
      pgn_bool_island[gd_num][pgnwh_num] = check_for_neighbours(gd_vector,
                                                                outer);
    }
  }

  /* Count and print total number of islands and non-islands. */
  int num_islands = 0;
  int num_non_islands = 0;
  for (int i = 0; i < (int) pgn_bool_island.size(); i++) {
    for (int j = 0; j < (int) pgn_bool_island[i].size(); j++) {
      num_islands = pgn_bool_island[i][j] ? num_islands + 1 : num_islands;
      num_non_islands = !pgn_bool_island[i][j] ? num_non_islands + 1
                                                 : num_non_islands; 
    }
  }
  std::cout << "Number of islands: " << num_islands << std::endl;
  std::cout << "Number of non-islands: " << num_non_islands << std::endl;

  return pgn_bool_island;
}

void insert_pll_into_graph(Polyline &pll, 
                           Graph &graph,
                           Point_vertex_map &pv_map) 
{
  vertex_descriptor u, v;
  for (int i = 0; i < (int) pll.size(); i++) {

    /* If the point is not yet in the graph. */
    if (pv_map.find(pll[i]) == pv_map.end()) {
      v = add_vertex(graph);
      pv_map[pll[i]] = v;
    } else {
      v = pv_map[pll[i]];
    }

    /* Associate the point to the vertex. */
    graph[v] = pll[i];  
    if (i != 0) {
      add_edge(u, v, graph);
    }
    u = v;
  }
}

/**** 3. Functions to create graph and split graph into unique polylines. ****/

Graph create_pll_graph(std::vector<GeoDiv> gd_vector)
{
  Graph graph;
  Point_vertex_map pv_map;
  for (GeoDiv gd : gd_vector) {
    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      Polyline pll_outer;
      Polygon outer = pgnwh.outer_boundary();
      for (Point pt_outer : outer) {
        pll_outer.push_back(pt_outer);
      }

      insert_pll_into_graph(pll_outer, graph, pv_map);

      std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
      for (Polygon hole : holes_v) {
        Polyline pll_hole;
        for (Point pt_hole : hole) {
          pll_hole.push_back(pt_hole);
        }

        insert_pll_into_graph(pll_hole, graph, pv_map);
      }
    }
  }
  return graph;
}

/* Struct template using Graph internally. */
template <typename Graph>
struct Polyline_visitor
{
  std::list<Polyline>& polyline_list;
  const Graph& points_pmap;

  /* Member initializer list. */
  Polyline_visitor(std::list<Polyline>& lines,
                   const Graph& points_property_map)
    : polyline_list(lines), points_pmap(points_property_map)
  {}

  void start_new_polyline() {
    Polyline V;
    polyline_list.push_back(V);
  }

  void add_node(typename boost::graph_traits<Graph>::vertex_descriptor vd) {
    Polyline& polyline = polyline_list.back();
    polyline.push_back(points_pmap[vd]);
  }

  void end_polyline() {}
};

void print_pll(PLL pll) {
  std::cout << pll.get_gd();
  std::cout << " | " << pll.get_pgnwh();
  std::cout << " | " << pll.get_pos();
  std::cout << " | " << pll.get_bool_hole();
  std::cout << " | " << pll.get_v1();
  std::cout << " | " << pll.get_vl();
  std::cout << " | " << pll.get_v2() << std::endl; 
}

/**
 * Track the matching progress of its polylines to its original gds and pgnwhs.
 */
class PgnProg {
  private:
    std::vector<Point> endpts;
    int island = false;
    int num_holes = 0;
    bool endpts_prog = false;
    bool holes_prog = false;

  public:
    void add_endpts(Point endpt1, Point endpt2) {
      endpts.push_back(endpt1);
      endpts.push_back(endpt2);
    }

    std::vector<Point> get_endpts() {
      return endpts;
    }

    void add_holes(int num_holes_) {
      num_holes += num_holes_;
    }

    void rem_hole() {
      num_holes--;
    }

    void update_prog() {
      if (num_holes == 0) {
        holes_prog = true;
      }

      std::sort(endpts.begin(), endpts.end());
      int pairs = 0;
      for (int i = 0; i < (int) endpts.size() - 1; i += 2) {
        if (endpts[i] == endpts[i + 1]) {
          pairs++;
        }
      }
      if (pairs * 2 == (int) endpts.size()) {
        endpts_prog = true;
      }
    }

    /** 
     * lone_pgn includes islands and polygons enclosed by a hole.
     * e.g. Brussels
     */
    void set_lone_pgn() {
      island = num_holes == 0 ? true : false;
    }

    bool check_prog() {
      return island || (endpts_prog && holes_prog);
    }
};

/**
 * Create a map where each pgnwh has an associated PgnProg object to track
 * the matching progress of its polylines to its original gds and pgnwhs.
 */
std::map<int, std::map<int, PgnProg>> create_pgnprog_map(
                                      std::vector<GeoDiv> gd_vector) 
{
  std::map<int, std::map<int, PgnProg>> pgnprog_map;
  for (int i = 0; i < (int) gd_vector.size(); i++) {
    for (int j = 0; j < (int) gd_vector[i].polygons_with_holes().size(); j++) {
      Polygon_with_holes pgnwh = gd_vector[i].polygons_with_holes()[j];
      std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
      pgnprog_map[i][j].add_holes(holes_v.size());
    }
  }
  return pgnprog_map;
}

int check_if_pll_on_pgn_boundary(PLL pll, Polygon pgn)
{
  /* Iterate through all vertices (points) inside the polyline. */
  int num_v_on_boundary = 0;
  for (int i = 0; i < (int) pll.get_pll().size(); i++) {
    Point pt = pll.get_pll()[i];
    bool v_on_boundary = CGAL::bounded_side_2(
                         pgn.begin(), pgn.end(), pt) == CGAL::ON_BOUNDARY;

    /* If the vertex is on the pgn's boundary, increment. */
    if (v_on_boundary) num_v_on_boundary++; 

    /** 
     * If the number of iterations has reached 3, break. This means that the
     * number of vertices on the pgn's boundary may or may not have reached 3.
     */
    if (i == 2) break;
  }

  return num_v_on_boundary;
}

void check_pgnprog_map(PLL pll,
    std::map<int, std::map<int, PgnProg>> &pgnprog_map,
    std::vector<int> &plls_match,
    bool pgnwh_is_island) {

  int pos = pll.get_pos();
  bool is_hole = pll.get_bool_hole();

  if (pll.get_v1() == pll.get_vl() && is_hole) {
    pgnprog_map[pll.get_gd()][pll.get_pgnwh()].rem_hole();
    plls_match[pos] += 1;
  } else if (pll.get_v1() == pll.get_vl() && pgnwh_is_island) {
    pgnprog_map[pll.get_gd()][pll.get_pgnwh()].set_lone_pgn();
    plls_match[pos] += 2;
  } else if (pll.get_v1() == pll.get_vl() && !is_hole) {
    pgnprog_map[pll.get_gd()][pll.get_pgnwh()].set_lone_pgn();
    plls_match[pos] += 1;
  } else {
    pgnprog_map[pll.get_gd()][pll.get_pgnwh()].add_endpts(pll.get_v1(), pll.get_vl());
    pgnprog_map[pll.get_gd()][pll.get_pgnwh()].update_prog();
    plls_match[pos] += 1;
  }
}

void check_prog_store_pll(int num_v_on_boundary,
                          PLL pll,
                          Polygon pgn,
                          std::map<int, std::vector<PLL>> &plls_by_pos,
                          std::map<int, std::map<int, PgnProg>> &pgnprog_map,
                          std::vector<int> &plls_match,
                          bool pgnwh_is_island)
{
  /**
   * If the number of vertices on the pgn's boundary reached 3:
   * a) Check pgnprog for this pgn taking into account this pll.
   * b) Store this pll according to its original map_state position
   *   along with its associated GeoDiv and Polygon_with_hole.
   */
  if (num_v_on_boundary >= 3) {
    check_pgnprog_map(pll, pgnprog_map, plls_match, pgnwh_is_island);
    plls_by_pos[pll.get_pos()].push_back(pll);

    /**
     * If the number of vertices on the pgn's boundary only reached 2, iterate
     * through all vertices of the pgn and determine if the 2 common vertices
     * in the pll and pgn are consecutive, checking both the clockwise and
     * counter-clockwise orientation of the pll and pgn.
     */
  } else if (num_v_on_boundary == 2) {
    for (int i = 0; i < (int) pgn.size() - 1; i++) {
      bool direction_1 = pgn[i] == pll.get_v1() && pgn[i + 1] == pll.get_v2();
      bool direction_2 = pgn[i + 1] == pll.get_v1() && pgn[i] == pll.get_v2();

      /** If the 2 common vertices are consecutive:
       * a) Check pgnprog for this pgn taking into account this pll.
       * b) Store this pll according to its original map_state position
       *    along with its associated GeoDiv and Polygon_with_hole.
       * c) Break because we no longer need to check all other vertices.
       */
      if (direction_1 || direction_2) {
        check_pgnprog_map(pll, pgnprog_map, plls_match, pgnwh_is_island);
        plls_by_pos[pll.get_pos()].push_back(pll);
        break;
      }
    }
  }
}

/** 6. Match and store polylines according to their original map_state      **/
/**    positions along with their associated GeoDiv and Polygon_with_hole.  **/

std::map<int, std::vector<PLL>> store_by_pos(
    std::vector<Polyline> &ct_polylines, 
    std::vector<GeoDiv> gd_vector,
    std::vector<std::vector<bool>> gd_pgnwh_island_bool)
{
  /**
   * Create map to match and store polylines according to their original
   * map_state positions along with their associated GeoDiv and
   * Polygon_with_hole.
   */
  std::map<int, std::vector<PLL>> plls_by_pos; 

  /**
   * Create a map where each pgnwh has an associated PgnProg object to track
   * the matching progress of its polylines to its original gds and pgnwhs.
   */
  std::map<int, std::map<int, PgnProg>> pgnprog_map = create_pgnprog_map(gd_vector);

  /**
   * Create a vector of polylines and an associated value x,
   * where x=1 or x=2, and x represents the number of polygons
   * it has been matched to.
   */
  std::vector<int> plls_match(ct_polylines.size(), 0);

  std::cout << "Matching polylines to polygons..." << std::endl;

  /* Iterate through all gds and pgnwhs and then through all polylines. */
  for (int gd_num = 0; gd_num < (int) gd_vector.size(); gd_num++) {
    std::cout << "gd: " << gd_num << std::endl;
    for (int pgnwh_num = 0;
         pgnwh_num < (int) gd_vector[gd_num].polygons_with_holes().size();
         pgnwh_num++) {

      for (int pos = 0; pos < (int) ct_polylines.size(); pos++) {
        Polyline polyl = ct_polylines[pos];

        /* If the pgnwh has already found/matched all its polylines, break */
        if (pgnprog_map[gd_num][pgnwh_num].check_prog()) break;

        /* If the polyl has been fully matched, continue */
        if (plls_match[pos] == 2) continue;

        Polygon_with_holes pgnwh = gd_vector[gd_num]
                                   .polygons_with_holes()[pgnwh_num];

        bool pgnwh_is_island = gd_pgnwh_island_bool[gd_num][pgnwh_num];

        /**
         * Create a PLL object with the current polyline's position,
         * Polyline object, gd_num, pgnwh_num, and whether the polyl
         * is a hole or not.
         */
        PLL pll_outer(pos, polyl, gd_num, pgnwh_num, false);
        Polygon outer_pgn = pgnwh.outer_boundary();

        /**
         * Check if pll_outer is on the boundary of outer_pgn and return 
         * the number of pll vertices on the boundary of the polygon.
         */
        int num_v_pll_outer = check_if_pll_on_pgn_boundary(pll_outer,
                                                           outer_pgn);

        /**
         * Depending on num_v_pll_outer, call check_pgnprog_map and plls_by_pos.
         */
        check_prog_store_pll(num_v_pll_outer,
                             pll_outer,
                             outer_pgn,
                             plls_by_pos,
                             pgnprog_map,
                             plls_match,
                             pgnwh_is_island);

        std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
        for (Polygon hole : holes_v) {
          PLL pll_hole(pos, polyl, gd_num, pgnwh_num, true);

          /* Check if pll_hole is on the boundary of hole. (for each hole) */
          int num_v_pll_hole = check_if_pll_on_pgn_boundary(pll_hole, hole);

          /**
           * Depending on num_v_pll_hole, call check_pgnprog_map and plls_by_pos.
           */
          check_prog_store_pll(num_v_pll_hole,
                               pll_hole,
                               hole,
                               plls_by_pos,
                               pgnprog_map,
                               plls_match,
                               pgnwh_is_island);
        }
      }
    }
  }
  std::cout << std::endl;

  return plls_by_pos;
}

std::map<int, std::map<int, std::vector<PLL>>> store_by_gd_pgnwh(std::vector<GeoDiv> container, 
    CT &ct, std::map<int, 
    std::vector<PLL>> &plls_by_pos) {
  std::map<int, std::map<int, std::vector<PLL>>> plls_by_gd_pgnwh;
  for (int gd_num = 0; gd_num < (int) container.size(); gd_num++) {
    for (int pgnwh_num = 0; pgnwh_num < (int) container[gd_num].polygons_with_holes().size(); pgnwh_num++) {
      int cit_num = 0;
      for (auto cit = ct.constraints_begin(); cit != ct.constraints_end(); cit++) {
        for (PLL pll : plls_by_pos[cit_num]) {
          if (pll.get_gd() == gd_num && pll.get_pgnwh() == pgnwh_num) {
            Polyline pll_ct;
            for (auto vit = ct.points_in_constraint_begin(*cit); vit != ct.points_in_constraint_end(*cit); vit++)
              pll_ct.push_back(*vit);
            PLL pll_new(pll.get_pos(), pll_ct, pll.get_gd(), pll.get_pgnwh(), pll.get_bool_hole());
            plls_by_gd_pgnwh[gd_num][pgnwh_num].push_back(pll_new);
            break;
          }
        }
        cit_num++;
      }
    }
  }
  return plls_by_gd_pgnwh;
}

void label_holes_correctly(std::vector<GeoDiv> container,
    std::map<int, std::map<int, std::vector<PLL>>> &plls_by_gd_pgnwh) {
  for (auto [gd_num, m] : plls_by_gd_pgnwh) {
    for (auto [pgnwh_num, v] : m) {
      for (PLL pll : v) {
        int gd_num2 = 0;
        for (GeoDiv gd3 : container) {
          int pgnwh_num2 = 0;
          for (Polygon_with_holes pgnwh3 : gd3.polygons_with_holes()) {
            std::vector<Polygon> holes_v(pgnwh3.holes_begin(), pgnwh3.holes_end());
            if (gd_num2 == gd_num && pgnwh_num2 == pgnwh_num && holes_v.empty()) {
              for (int i = 0; i < (int) plls_by_gd_pgnwh[gd_num][pgnwh_num].size(); i++)
                pll.set_bool_hole(false);
            } else {
              break;
            }
            pgnwh_num2++;
          }
          gd_num2++;
        }
      }
    }
  }
}

void set_visited_vals(std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, bool>>> &visited,
    std::map<int, std::map<int, std::vector<PLL>>> &plls_by_gd_pgnwh) {
  for (auto [gd_num, m] : plls_by_gd_pgnwh) {
    for (auto [pgnwh_num, pll_v] : m) {
      int i = 0;
      for (PLL pll : pll_v) {
        visited[gd_num][pgnwh_num][pll.get_pos()] = false; // Set all visited to false
        i++;
      }
      // Sort each vector<PLL> so that all holes are in front
      std::sort(plls_by_gd_pgnwh[gd_num][pgnwh_num].begin(), plls_by_gd_pgnwh[gd_num][pgnwh_num].end(), 
          [](PLL pll1, PLL pll2) {
          return pll1.get_bool_hole() > pll2.get_bool_hole();
          });
    }
  }
  // Print out new sequence
  /*
     for (auto [gd_num, m] : plls_by_gd_pgnwh)
     for (auto [pgnwh_num, pll_v] : m)
     for (PLL pll : pll_v)
     print_pll(pll);
     std::cout << std::endl;
     */
}

void assemble_pll_to_pgn(std::map<int, std::map<int, std::vector<PLL>>> &plls_by_gd_pgnwh, 
    std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, bool>>> &visited, 
    std::vector<GeoDiv> &container_final,
    std::vector<GeoDiv> container_org)
{
  std::cout << "Assembling polylines into polygons..." << std::endl;
  for (auto [gd_num, m] : plls_by_gd_pgnwh) {
    GeoDiv gd_final(container_org[gd_num].id());
    std::vector<Polygon> holes_v;
    for (auto [pgnwh_num, pll_v] : m) {

      for (PLL pll : pll_v) {
        Polygon outer; // This will only be for islands/holes anyway
        if (visited[gd_num][pgnwh_num][pll.get_pos()]) continue;

        // std::cout << pll.get_gd() << " " << pll.get_pgnwh() << " " << pll.get_pos() << " " << pll.get_bool_hole() << std::endl;

        // if it is a single polyline (e.g. island)
        if (pll.get_v1() == pll.get_vl() && !visited[gd_num][pgnwh_num][pll.get_pos()]) {

          visited[gd_num][pgnwh_num][pll.get_pos()] = true;

          for (Point pt : pll.get_pll())
            outer.push_back(pt);

          if (!pll.get_bool_hole() && holes_v.empty()) {
            Polygon_with_holes pgnwh(outer);
            gd_final.push_back(pgnwh);
          } else if (!pll.get_bool_hole() && !holes_v.empty()) { // if there is a hole within
            // Check if hole's middle vertex is inside boundary
            bool holes_inside = true; 

            for (Polygon hole : holes_v)
              if (CGAL::bounded_side_2(outer.begin(), outer.end(), hole[hole.size() / 2]) != CGAL::ON_BOUNDED_SIDE)
                holes_inside = false;

            if (holes_inside) {
              Polygon_with_holes pgnwh_final_2(outer, holes_v.begin(), holes_v.end());
              gd_final.push_back(pgnwh_final_2);
              holes_v.clear();
            }
          } else {
            holes_v.push_back(outer);
          }
        } else {
          std::deque<PLL> deq;
          deq.push_back(pll);
          visited[gd_num][pgnwh_num][pll.get_pos()] = true;

          // Connect together all polylines belonging to a Polygon_with_holes
          while (1) {
            bool found = false;

            for (PLL pll2 : plls_by_gd_pgnwh[pll.get_gd()][pll.get_pgnwh()]) {
              if (visited[gd_num][pgnwh_num][pll2.get_pos()] == true) continue;

              if (deq.front().get_v1()[0] == pll2.get_vl()[0] &&
                  deq.front().get_v1()[1] == pll2.get_vl()[1]) {
                deq.push_front(pll2);
                visited[gd_num][pgnwh_num][pll2.get_pos()] = true;
                found = true;
              } else if (deq.back().get_vl()[0] == pll2.get_v1()[0] &&
                  deq.back().get_vl()[1] == pll2.get_v1()[1]) {
                deq.push_back(pll2);
                visited[gd_num][pgnwh_num][pll2.get_pos()] = true;
                found = true;
              } else if (deq.front().get_v1()[0] == pll2.get_v1()[0] &&
                  deq.front().get_v1()[1] == pll2.get_v1()[1]) {
                Polyline polyl_new = pll2.get_pll();
                std::reverse(polyl_new.begin(), polyl_new.end());
                pll2.set_pll(polyl_new);
                deq.push_front(pll2);
                visited[gd_num][pgnwh_num][pll2.get_pos()] = true;
                found = true;
              } else if (deq.back().get_vl()[0] == pll2.get_vl()[0] &&
                  deq.back().get_vl()[1] == pll2.get_vl()[1]) {
                Polyline polyl_new = pll2.get_pll();
                std::reverse(polyl_new.begin(), polyl_new.end());
                pll2.set_pll(polyl_new);
                deq.push_back(pll2);
                visited[gd_num][pgnwh_num][pll2.get_pos()] = true;
                found = true;
              }
            }
            if (!found) break;
          }
          Polygon outer_2;

          for (PLL pll_deq : deq)
            for (Point pt : pll_deq.get_pll())
              outer_2.push_back(pt);

          if (pll.get_bool_hole()) {
            holes_v.push_back(outer_2);
            continue;
          }

          if (holes_v.empty()) {
            Polygon_with_holes pgnwh_final_2(outer_2);
            gd_final.push_back(pgnwh_final_2);
          } else {
            // Check if hole's middle vertex is inside boundary
            bool holes_inside = true; 
            for (Polygon hole : holes_v)
              if (CGAL::bounded_side_2(outer_2.begin(), outer_2.end(), hole[hole.size() / 2]) != CGAL::ON_BOUNDED_SIDE)
                holes_inside = false;
            if (holes_inside) {
              Polygon_with_holes pgnwh_final_2(outer_2, holes_v.begin(), holes_v.end());
              gd_final.push_back(pgnwh_final_2);
              holes_v.clear();
            }
          }
          // std::cout << std::endl;
        }

      }
    }
    container_final.push_back(gd_final);
  }
}

void print_num_pts(std::vector<GeoDiv> container) {
  int num_pts = 0;
  for (GeoDiv gd : container) {
    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      Polygon outer = pgnwh.outer_boundary();
      num_pts += outer.size();

      std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
      for (Polygon hole : holes_v)
        num_pts += hole.size();
    }
  }
  std::cout << num_pts << std::endl;
}

void remove_first_point_as_last_point(std::vector<GeoDiv> &container) {
  for (GeoDiv &gd : container) {
    for (Polygon_with_holes &pgnwh : *gd.ref_to_polygons_with_holes()) {

      Polygon *outer = &pgnwh.outer_boundary();
      auto outer_it = outer->end();
      outer_it--;
      outer->erase(outer_it);

      std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
      for (auto hole = pgnwh.holes_begin(); hole != pgnwh.holes_end(); hole++) {
        auto hole_it = hole->end();
        hole_it--;
        hole->erase(hole_it);
      }
    }
  }
}

void simplify_map(MapState *map_state)
{
  /********************************* Steps: **********************************/
  /* 1. Function to repeat the first point as the last point.                */
  /* 2. Function to return a vector<vector<bool>> of pgns IDed as islands.   */
  /* 3. Functions to create graph and split graph into unique polylines.     */ 
  /* 4. Store polylines in a CT object from polyline_list.                   */
  /* 5. Store polylines in a vector in the order of CT polylines.            */
  /* 6. Match and store polylines according to their original map_state      */
  /*    positions along with their associated GeoDiv and Polygon_with_hole.  */
  /* 7. Simplify polylines.                                                  */
  // 8. Store polylines by GeoDivs and Polygon_with_holes with their associated positions
  // 9. Set visited values
  // 10. Assemble polylines into polygons
  // 11. Remove first point as last point by reference 

  std::vector<GeoDiv> container = map_state->geo_divs();

  /* 1. Function to repeat the first point as the last point.                */
  repeat_first_point_as_last_point(container);

  /* Start timer for step 2. */
  const auto start_s2 = std::chrono::system_clock::now();

  /* 2. Function to return a vector<vector<bool>> of pgns IDed as islands    */
  std::vector<std::vector<bool>> pgn_bool_island = 
    get_pgn_bool_island(container);

  /* Get and print step 2's elapsed time */
  const std::chrono::duration<double, std::milli> dur_s2 =
        std::chrono::system_clock::now() - start_s2;
  std::cout << "get_gd_pgnwh_island_bool() time elapsed: ";
  std::cout << dur_s2.count() << " ms (";
  std::cout << dur_s2.count() / 1000 << " s)" << std::endl;
  std::cout << std::endl;

  /* Start timer for remaining steps. */
  const auto start_s311 = std::chrono::system_clock::now();

  /* 3. Functions to create graph and split graph into unique polylines.     */
  Graph graph = create_pll_graph(container);
  std::list<Polyline> polyline_list;
  Polyline_visitor<Graph> polyline_visitor(polyline_list, graph);
  CGAL::split_graph_into_polylines(graph, polyline_visitor);

  /* 4. Store polylines from polyline_list in a CT object.                   */
  CT ct; 
  for (Polyline polyline : polyline_list) {
    ct.insert_constraint(polyline.begin(), polyline.end());
  }

  /* 5. Store polylines in a vector in the order of CT polylines.            */
  /**
   * This is so that we no longer need to use the CT data structure when
   * carrying out operations, which has limited methods available to use.
   * Instead, we can use the vector data structure and its methods.
   */
  std::vector<Polyline> ct_polylines;
  for (auto cit = ct.constraints_begin(); cit != ct.constraints_end()
                                        ; cit++) {
    Polyline polyl;
    for (auto vit = ct.points_in_constraint_begin(*cit)
         ; vit != ct.points_in_constraint_end(*cit); vit++) {
      polyl.push_back(*vit);
    }
    ct_polylines.push_back(polyl);
  }

  /* 6. Match and store polylines according to their original map_state      */
  /*    positions along with their associated GeoDiv and Polygon_with_hole.  */
  std::map<int, std::vector<PLL>> 
    plls_by_pos = store_by_pos(ct_polylines, container, pgn_bool_island);

  /* 7. Simplify polylines.                                                  */
  PS::simplify(ct, Cost(), Stop(0.2));

  // 8. Store polylines by GeoDivs and Polygon_with_holes with their associated positions
  std::map<int, std::map<int, std::vector<PLL>>> plls_by_gd_pgnwh = store_by_gd_pgnwh(container, ct, plls_by_pos);

  // 9. Set visited values
  std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, bool>>> visited;
  set_visited_vals(visited, plls_by_gd_pgnwh);

  // 10. Assemble polylines into polygons
  std::vector<GeoDiv> container_simp;
  assemble_pll_to_pgn(plls_by_gd_pgnwh, visited, container_simp, container);

  // Print number of points before simplifying (container)
  std::cout << "Number of vertices before simplifying (container): ";
  print_num_pts(container);

  // Print number of points after simplifying
  std::cout << "Number of vertices after simplifying: ";
  print_num_pts(container_simp);

  // 11. Remove first point as last point by reference 
  remove_first_point_as_last_point(container_simp);

  map_state->set_geo_divs(container_simp);

  const std::chrono::duration<double, std::milli> dur_s311 = std::chrono::system_clock::now() - start_s311;
  std::cout << "Remaining simplification steps time elapsed: " << dur_s311.count() << " ms (";
  std::cout << dur_s311.count() / 1000 << " s)" << std::endl;
  std::cout << std::endl;

  double dur_total = dur_s2.count() + dur_s311.count();
  std::cout << "Total simplification time elapsed: " << dur_total << " ms (";
  std::cout << dur_total / 1000 << " s)" << std::endl;
  std::cout << std::endl;
}
