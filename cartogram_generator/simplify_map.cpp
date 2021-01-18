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

// Inserts a polyline into the graph
void insert(const std::vector<Point>& poly, 
            Graph& graph,
            Point_vertex_map& pvmap) {
  vertex_descriptor u, v;
  for (int i = 0; i < (int) poly.size(); i++) {
    // check if the point is not yet in the graph
    if (pvmap.find(poly[i]) == pvmap.end()) {
      v = add_vertex(graph);
      pvmap[poly[i]] = v;
    } else {
      v = pvmap[poly[i]];
    }
    graph[v] = poly[i];  // associate the point to the vertex
    if (i != 0) {
      add_edge(u, v, graph);
    }
    u = v;
  }
}

template <typename Graph>
struct Polyline_visitor {
  std::list<Polyline>& polylines;
  const Graph& points_pmap;

  Polyline_visitor(std::list<Polyline>& lines, const Graph& points_property_map)
    : polylines(lines), points_pmap(points_property_map)
  {}

  void start_new_polyline() {
    Polyline V;
    polylines.push_back(V);
  }

  void add_node(typename boost::graph_traits<Graph>::vertex_descriptor vd) {
    Polyline& polyline = polylines.back();
    polyline.push_back(points_pmap[vd]);
  }

  void end_polyline() {}
};

Graph create_pll_graph(std::vector<GeoDiv> container) {
  Graph graph;
  Point_vertex_map pvmap;
  for (GeoDiv gd : container) {
    for (Polygon_with_holes pgnwh : gd.polygons_with_holes()) {
      Polyline pll_outer; // polyline for outer boundary
      Polygon outer = pgnwh.outer_boundary();
      for (Point pt_outer : outer)
        pll_outer.push_back(pt_outer);
      insert(pll_outer, graph, pvmap); // insert pll_outer into graph

      std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
      for (Polygon hole : holes_v) {
        Polyline pll_hole; // polyline for hole
        for (Point pt_hole : hole)
          pll_hole.push_back(pt_hole);
        insert(pll_hole, graph, pvmap); // insert pll_hole into graph
      }
    }
  }
  return graph;
}

void print_pll(PLL pll) {
  std::cout << pll.get_gd();
  std::cout << " | " << pll.get_pgnwh();
  std::cout << " | " << pll.get_pos();
  std::cout << " | " << pll.get_bool_hole();
  std::cout << " | " << pll.get_v1();
  std::cout << " | " << pll.get_vl();
  std::cout << " | " << pll.get_v2() << std::endl; 
}

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
      for (int i = 0; i < (int) endpts.size() - 1; i += 2)
        if (endpts[i] == endpts[i + 1])
          pairs++;
      if (pairs * 2 == (int) endpts.size())
        endpts_prog = true;
    }

    // lone_pgn includes islands and polygons enclosed by a hole (e.g. Brussels)
    void set_lone_pgn() {
      island = num_holes == 0 ? true : false;
    }

    bool check_prog() {
      return island || (endpts_prog && holes_prog);
    }
};

std::map<int, std::map<int, PgnProg>> create_pgn_prog(std::vector<GeoDiv> container_dens) {
  std::map<int, std::map<int, PgnProg>> pgn_prog;
  for (int i = 0; i < (int) container_dens.size(); i++) {
    for (int j = 0; j < (int) container_dens[i].polygons_with_holes().size(); j++) {
      Polygon_with_holes pgnwh = container_dens[i].polygons_with_holes()[j];
      std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
      pgn_prog[i][j].add_holes(holes_v.size());
    }
  }
  return pgn_prog;
}

void check_pgn_prog(PLL pll,
    std::map<int, std::map<int, PgnProg>> &pgn_prog,
    std::vector<int> &matched_pgns) {

  int pos = pll.get_pos();
  bool is_hole = pll.get_bool_hole();

  if (pll.get_v1() == pll.get_vl() && is_hole) {
    pgn_prog[pll.get_gd()][pll.get_pgnwh()].rem_hole();
    matched_pgns[pos] += 1;
  } else if (pll.get_v1() == pll.get_vl() && !is_hole) {
    pgn_prog[pll.get_gd()][pll.get_pgnwh()].set_lone_pgn();
    matched_pgns[pos] += 1;
  } else {
    pgn_prog[pll.get_gd()][pll.get_pgnwh()].add_endpts(pll.get_v1(), pll.get_vl());
    pgn_prog[pll.get_gd()][pll.get_pgnwh()].update_prog();
    matched_pgns[pos] += 1;
  }
}

void check_if_pll_on_pgn_boundary(PLL pll,
    Polygon pgn,
    std::map<int, std::vector<PLL>> &pll_cntr_by_pos,
    std::map<int, std::map<int, PgnProg>> &pgn_prog,
    std::vector<int> &matched_pgns) {

  // need >=3 vertices in pll to be on pgn's boundary to count as part of pgn
  int num_v_on_outer = 0;
  for (int i = 0; i < (int) pll.get_pll().size(); i++) {
    Point pt = pll.get_pll()[i];
    bool v_on_outer = CGAL::bounded_side_2(pgn.begin(), pgn.end(), pt) == CGAL::ON_BOUNDARY;
    if (v_on_outer) num_v_on_outer++; 
    if (i == 2) break;
  }

  int pos = pll.get_pos();

  if (num_v_on_outer >= 3) {
    // print_pll(pll);

    check_pgn_prog(pll, pgn_prog, matched_pgns);
    pll_cntr_by_pos[pos].push_back(pll);

  } else if (num_v_on_outer == 2) {
    // Consecutive pll check: accounts for clockwise and counter-clockwise orientation of pgn
    for (int i = 0; i < (int) pgn.size() - 1; i++) {
      bool direction_1 = pgn[i] == pll.get_v1() && pgn[i + 1] == pll.get_v2();
      bool direction_2 = pgn[i + 1] == pll.get_v1() && pgn[i] == pll.get_v2();
      if (direction_1 || direction_2) {
        // print_pll(pll);

        check_pgn_prog(pll, pgn_prog, matched_pgns);
        pll_cntr_by_pos[pos].push_back(pll);
        break;
      }
    }
  }
}

bool are_equal(Polyline polyl1, Polyline polyl2) {
  if (polyl1.size() == 0 || polyl2.size() == 0) return false;

  bool same_vertices = polyl1[0] == polyl2[0]
    && (polyl1[1] == polyl2[1] || polyl1[1] == polyl2[polyl2.size() - 2])
    && (polyl1[2] == polyl2[2] || polyl1[2] == polyl2[polyl2.size() - 3])
    && polyl1[polyl1.size() - 1] == polyl2[polyl2.size() - 1];

  return same_vertices;
}

bool is_duplicate(std::vector<GeoDiv> container_dens,
    std::map<int, std::map<int, std::pair<Polyline, bool>>> is_island,
    Polyline polyl) {

  for (int gd_num = 0; gd_num < (int) container_dens.size(); gd_num++) {
    for (int pgnwh_num = 0; pgnwh_num < (int) container_dens[gd_num].polygons_with_holes().size(); pgnwh_num++) {
      if (are_equal(is_island[gd_num][pgnwh_num].first, polyl))
        return true;
    }
  }
  return false;
}

bool update_duplicated(std::map<int, std::map<int, std::pair<Polyline, bool>>> &is_island) {
}

void identify_islands(std::vector<GeoDiv> container_dens,
    std::vector<Polyline> ct_polylines,
    std::map<int, std::map<int, std::pair<Polyline, bool>>> &is_island) {

  for (int gd_num = 0; gd_num < (int) container_dens.size(); gd_num++)
    for (int pgnwh_num = 0; pgnwh_num < (int) container_dens[gd_num].polygons_with_holes().size(); pgnwh_num++)
      is_island[gd_num][pgnwh_num] = {{}, false};

  for (int gd_num = 0; gd_num < (int) container_dens.size(); gd_num++) {
    std::cout << "gd: " << gd_num << std::endl;
    for (int pgnwh_num = 0; pgnwh_num < (int) container_dens[gd_num].polygons_with_holes().size(); pgnwh_num++) {

      bool match = false;
      Polygon_with_holes pgnwh = container_dens[gd_num].polygons_with_holes()[pgnwh_num];

      // outer
      Polygon outer = pgnwh.outer_boundary();
      Polyline polyl_outer;
      for (Point pt_outer : outer) polyl_outer.push_back(pt_outer);
      for (Polyline polyl : ct_polylines) {
        if (are_equal(polyl_outer, polyl)) {
          match = true;
          if (is_duplicate(container_dens, is_island, polyl)) {
            std::cout << "duplicate" << std::endl;
            is_island[gd_num][pgnwh_num] = {polyl, false};
          } else {
            is_island[gd_num][pgnwh_num] = {polyl, true};
          }
          std::cout << "match outer" << std::endl; 
          break;
        }
      }

      // hole 
      if (!match) {
        std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
        for (Polygon hole : holes_v) {
          Polyline polyl_hole;
          for (Point pt_hole: hole) polyl_hole.push_back(pt_hole);
          for (Polyline polyl : ct_polylines) {
            if (are_equal(polyl_hole, polyl)) {
              match = true;
              if (is_duplicate(container_dens, is_island, polyl)) {
                std::cout << "duplicate" << std::endl;
                is_island[gd_num][pgnwh_num] = {polyl, false};
              } else {
                is_island[gd_num][pgnwh_num] = {polyl, true};
              }
              std::cout << "match hole" << std::endl; 
              break;
            }
          }
          if (match) break;
        }
      }
    }
  }
}

std::map<int, std::vector<PLL>> store_by_pos(std::vector<Polyline> &ct_polylines, 
    std::vector<GeoDiv> container_dens) {

  // Create map to check if geo_divs and pgnwhs have been properly matched
  std::map<int, std::vector<PLL>> pll_cntr_by_pos; 

  // Create map to track polygon matching progress
  std::map<int, std::map<int, PgnProg>> pgn_prog = create_pgn_prog(container_dens);

  // Create vector of visited polylines
  std::vector<int> matched_pgns(ct_polylines.size(), 0);

  for (int gd_num = 0; gd_num < (int) container_dens.size(); gd_num++) {
    std::cout << "gd: " << gd_num << std::endl;
    for (int pgnwh_num = 0; pgnwh_num < (int) container_dens[gd_num].polygons_with_holes().size(); pgnwh_num++) {

      for (int pos = 0; pos < (int) ct_polylines.size(); pos++) {
        Polyline polyl = ct_polylines[pos];

        // If the pgnwh has already found all its polylines, break 
        if (pgn_prog[gd_num][pgnwh_num].check_prog()) break;

        // If the polyl has already been matched, continue
        if (matched_pgns[pos] == 2) continue;

        Polygon_with_holes pgnwh = container_dens[gd_num].polygons_with_holes()[pgnwh_num];
        PLL pll_outer(pos, polyl, gd_num, pgnwh_num, false);
        Polygon outer = pgnwh.outer_boundary();

        // Check outer polygon
        check_if_pll_on_pgn_boundary(pll_outer,
            outer,
            pll_cntr_by_pos,
            pgn_prog,
            matched_pgns);

        std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
        for (Polygon hole : holes_v) {
          PLL pll_hole(pos, polyl, gd_num, pgnwh_num, true);

          // Check holes
          check_if_pll_on_pgn_boundary(pll_hole,
              hole,
              pll_cntr_by_pos,
              pgn_prog,
              matched_pgns);
        }
      }
    }
  }
  std::cout << std::endl;

  return pll_cntr_by_pos;
}

std::map<int, std::map<int, std::vector<PLL>>> store_by_gd_pgnwh(std::vector<GeoDiv> container, 
    CT &ct, std::map<int, 
    std::vector<PLL>> &pll_cntr_by_pos) {
  std::map<int, std::map<int, std::vector<PLL>>> pll_cntr_by_gd_pgnwh;
  for (int gd_num = 0; gd_num < (int) container.size(); gd_num++) {
    for (int pgnwh_num = 0; pgnwh_num < (int) container[gd_num].polygons_with_holes().size(); pgnwh_num++) {
      int cit_num = 0;
      for (auto cit = ct.constraints_begin(); cit != ct.constraints_end(); cit++) {
        for (PLL pll : pll_cntr_by_pos[cit_num]) {
          if (pll.get_gd() == gd_num && pll.get_pgnwh() == pgnwh_num) {
            Polyline pll_ct;
            for (auto vit = ct.points_in_constraint_begin(*cit); vit != ct.points_in_constraint_end(*cit); vit++)
              pll_ct.push_back(*vit);
            PLL pll_new(pll.get_pos(), pll_ct, pll.get_gd(), pll.get_pgnwh(), pll.get_bool_hole());
            pll_cntr_by_gd_pgnwh[gd_num][pgnwh_num].push_back(pll_new);
            break;
          }
        }
        cit_num++;
      }
    }
  }
  return pll_cntr_by_gd_pgnwh;
}

void label_holes_correctly(std::vector<GeoDiv> container,
    std::map<int, std::map<int, std::vector<PLL>>> &pll_cntr_by_gd_pgnwh) {
  for (auto [gd_num, m] : pll_cntr_by_gd_pgnwh) {
    for (auto [pgnwh_num, v] : m) {
      for (PLL pll : v) {
        int gd_num2 = 0;
        for (GeoDiv gd3 : container) {
          int pgnwh_num2 = 0;
          for (Polygon_with_holes pgnwh3 : gd3.polygons_with_holes()) {
            std::vector<Polygon> holes_v(pgnwh3.holes_begin(), pgnwh3.holes_end());
            if (gd_num2 == gd_num && pgnwh_num2 == pgnwh_num && holes_v.empty()) {
              for (int i = 0; i < (int) pll_cntr_by_gd_pgnwh[gd_num][pgnwh_num].size(); i++)
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
    std::map<int, std::map<int, std::vector<PLL>>> &pll_cntr_by_gd_pgnwh) {
  for (auto [gd_num, m] : pll_cntr_by_gd_pgnwh) {
    for (auto [pgnwh_num, pll_v] : m) {
      int i = 0;
      for (PLL pll : pll_v) {
        visited[gd_num][pgnwh_num][pll.get_pos()] = false; // Set all visited to false
        i++;
      }
      // Sort each vector<PLL> so that all holes are in front
      std::sort(pll_cntr_by_gd_pgnwh[gd_num][pgnwh_num].begin(), pll_cntr_by_gd_pgnwh[gd_num][pgnwh_num].end(), 
          [](PLL pll1, PLL pll2) {
          return pll1.get_bool_hole() > pll2.get_bool_hole();
          });
    }
  }
  // Print out new sequence
  /*
     for (auto [gd_num, m] : pll_cntr_by_gd_pgnwh)
     for (auto [pgnwh_num, pll_v] : m)
     for (PLL pll : pll_v)
     print_pll(pll);
     std::cout << std::endl;
     */
}

void assemble_pll_to_pgn(std::map<int, std::map<int, std::vector<PLL>>> &pll_cntr_by_gd_pgnwh, 
    std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, bool>>> &visited, 
    std::vector<GeoDiv> &container_final,
    std::vector<GeoDiv> container_org) {
  for (auto [gd_num, m] : pll_cntr_by_gd_pgnwh) {
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

            for (PLL pll2 : pll_cntr_by_gd_pgnwh[pll.get_gd()][pll.get_pgnwh()]) {
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

void repeat_first_point_as_last_point(std::vector<GeoDiv> &container) {
  for (GeoDiv &gd : container) {
    for (Polygon_with_holes &pgnwh : *gd.ref_to_polygons_with_holes()) {
      Polygon_with_holes pgnwh_new;

      Polygon *outer = &pgnwh.outer_boundary();
      outer->push_back((*outer)[0]);

      std::vector<Polygon> holes_v(pgnwh.holes_begin(), pgnwh.holes_end());
      for (auto hole = pgnwh.holes_begin(); hole != pgnwh.holes_end(); hole++) {
        hole->push_back((*hole)[0]); 
      }
    }
  }
}

void remove_first_point_as_last_point(std::vector<GeoDiv> &container) {
  for (GeoDiv &gd : container) {
    for (Polygon_with_holes &pgnwh : *gd.ref_to_polygons_with_holes()) {
      Polygon_with_holes pgnwh_new;

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

void simplify_map(MapState *map_state) {
  // Steps:
  // 1. Repeat first point as last point by reference 
  // 2. Create graph and split graph into unique polylines
  // 3. Densify
  // 4. Store polylines from polyline_list in CT
  // 5. Store ct polylines (densified) with their associated polylines (non-densified)
  // 6. Store polylines by positions with their associated GeoDivs and Polygon_with_holes
  // 7. Simplify polylines
  // 8. Store polylines by GeoDivs and Polygon_with_holes with their associated positions
  // 9. Set visited values
  // 10. Assemble polylines into polygons
  // 11. Remove first point as last point by reference 

  const auto start = std::chrono::system_clock::now();

  std::vector<GeoDiv> container = map_state->geo_divs();

  // 1. Repeat first point as last point by reference 
  repeat_first_point_as_last_point(container);

  // Densify container
  std::vector<GeoDiv> container_dens = densify(container);

  // 2. Create graph and split graph into unique polylines
  Graph graph = create_pll_graph(container_dens);
  std::list<Polyline> polyline_list;
  Polyline_visitor<Graph> polyline_visitor(polyline_list, graph);
  CGAL::split_graph_into_polylines(graph, polyline_visitor);

  // 4. Store polylines from polyline_list in CT
  CT ct; 
  for (Polyline polyline : polyline_list)
    ct.insert_constraint(polyline.begin(), polyline.end());

  // Create vector to store polylines in CT order
  std::vector<Polyline> ct_polylines;
  for (auto cit = ct.constraints_begin(); cit != ct.constraints_end(); cit++) {
    Polyline polyl;
    for (auto vit = ct.points_in_constraint_begin(*cit); vit != ct.points_in_constraint_end(*cit); vit++)
      polyl.push_back(*vit);
    ct_polylines.push_back(polyl);
  }

  std::map<int, std::map<int, std::pair<Polyline, bool>>> is_island;
  identify_islands(container_dens, ct_polylines, is_island);
  const std::chrono::duration<double, std::milli> duration = std::chrono::system_clock::now() - start;
  std::cout << "simplify_map() time elapsed: " << duration.count() << "ms (";
  std::cout << duration.count() / 1000 << "s)" << std::endl;
  std::cout << std::endl;

  // 6. Store polylines by positions with their associated GeoDivs and Polygon_with_holes
  // std::cout << "Store polylines by positions with their associated GeoDivs and Polygon_with_holes" << std::endl;
  std::map<int, std::vector<PLL>> pll_cntr_by_pos = store_by_pos(ct_polylines, container_dens);

  // 7. Simplify polylines
  PS::simplify(ct, Cost(), Stop(0.2));

  // 8. Store polylines by GeoDivs and Polygon_with_holes with their associated positions
  // std::cout << "Store polylines by GeoDivs and Polygon_with_holes with their associated positions" << std::endl;
  std::map<int, std::map<int, std::vector<PLL>>> pll_cntr_by_gd_pgnwh = store_by_gd_pgnwh(container, ct, pll_cntr_by_pos);

  // 9. Set visited values
  std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, bool>>> visited;
  set_visited_vals(visited, pll_cntr_by_gd_pgnwh);

  // std::cout << "Assemble polylines into polygons" << std::endl;
  // 10. Assemble polylines into polygons
  std::vector<GeoDiv> container_simp;
  assemble_pll_to_pgn(pll_cntr_by_gd_pgnwh, visited, container_simp, container_dens);

  // Print number of points before simplifying (container)
  std::cout << "Number of vertices before simplifying (container): ";
  print_num_pts(container);

  // Print number of points before simplifying (container_dens)
  std::cout << "Number of vertices before simplifying (container_dens): ";
  print_num_pts(container_dens);

  // Print number of points after simplifying
  std::cout << "Number of vertices after simplifying: ";
  print_num_pts(container_simp);

  // 11. Remove first point as last point by reference 
  remove_first_point_as_last_point(container_simp);

  map_state->set_geo_divs(container_simp);

  const std::chrono::duration<double, std::milli> duration2 = std::chrono::system_clock::now() - start;
  std::cout << "simplify_map() time elapsed: " << duration2.count() << "ms (";
  std::cout << duration2.count() / 1000 << "s)" << std::endl;
}
