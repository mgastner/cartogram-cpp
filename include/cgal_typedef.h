#ifndef CGAL_TYPEDEF_H_
#define CGAL_TYPEDEF_H_

#include <CGAL/Barycentric_coordinates_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Quadtree.h>

typedef CGAL::Simple_cartesian<double> Scd;
typedef CGAL::Polygon_2<Scd> Polygon;
typedef CGAL::Polygon_with_holes_2<Scd> Polygon_with_holes;
typedef CGAL::Aff_transformation_2<Scd> Transformation;
typedef CGAL::Point_2<Scd> Point;
typedef CGAL::Line_2<Scd> Line;
typedef CGAL::Bbox_2 Bbox;
typedef CGAL::Segment_2<Scd> Segment;

// Minimum enclosing ellipses and circles
typedef CGAL::Min_ellipse_2_traits_2<Scd> Ellipse_traits;
typedef CGAL::Min_ellipse_2<Ellipse_traits> Min_ellipse;
typedef CGAL::Min_circle_2_traits_2<Scd> Circle_traits;
typedef CGAL::Min_circle_2<Circle_traits> Min_circle;

// Quadtree
typedef CGAL::Orthtree<CGAL::Orthtree_traits_2<Scd>, std::vector<Point>>
  Quadtree;

// Delaunay triangulation
typedef CGAL::Delaunay_triangulation_2<Scd> Delaunay;
typedef Delaunay::Line_face_circulator Line_face_circulator;
typedef Delaunay::Face_handle Face_handle;

// Polyline simplification
typedef CGAL::Polyline_simplification_2::Vertex_base_2<Scd> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<Scd> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef CGAL::
  Constrained_Delaunay_triangulation_2<Scd, TDS, CGAL::Exact_predicates_tag>
    CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT> CT;
typedef CT::Constraint_id Constraint_id;
typedef CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold Stop;
typedef CGAL::Polyline_simplification_2::Squared_distance_cost Cost;

#endif
