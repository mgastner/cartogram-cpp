#ifndef CGAL_TYPEDEF_HPP_
#define CGAL_TYPEDEF_HPP_

#include <CGAL/Barycentric_coordinates_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>

typedef CGAL::Simple_cartesian<double> Scd;
typedef CGAL::Polygon_2<Scd> Polygon;
typedef CGAL::Polygon_with_holes_2<Scd> Polygon_with_holes;
typedef CGAL::Aff_transformation_2<Scd> Transformation;
typedef CGAL::Point_2<Scd> Point;
typedef CGAL::Vector_2<Scd> Vector;
typedef CGAL::Line_2<Scd> Line;
typedef CGAL::Iso_rectangle_2<Scd> Bbox;
typedef CGAL::Segment_2<Scd> Segment;

#endif  // CGAL_TYPEDEF_HPP_
