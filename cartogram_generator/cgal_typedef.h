#ifndef CGAL_TYPEDEF_H_
#define CGAL_TYPEDEF_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <boost/graph/adjacency_list.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes;
typedef CGAL::Aff_transformation_2<K> Transformation;
typedef CGAL::Point_2<K> Point;

// Polyline
namespace PS = CGAL::Polyline_simplification_2;
typedef PS::Vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, CGAL::Exact_predicates_tag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT> CT;
typedef PS::Stop_below_count_ratio_threshold Stop;
typedef PS::Squared_distance_cost Cost;

// Split graph into polylines
typedef boost::adjacency_list<boost::vecS, boost::setS, boost::undirectedS, Point> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
typedef std::map<Point, vertex_descriptor> Point_vertex_map;
typedef std::vector<Point> Polyline;

#endif
