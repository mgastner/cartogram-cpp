#ifndef CGAL_TYPEDEF_H_
#define CGAL_TYPEDEF_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_with_holes_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;
typedef CGAL::Polygon_2<Epick> Polygon;
typedef CGAL::Polygon_with_holes_2<Epick> Polygon_with_holes;
typedef CGAL::Aff_transformation_2<Epick> Transformation;
typedef CGAL::Point_2<Epick> Point;

#endif
