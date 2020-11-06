#ifndef CGAL_TYPEDEF_H_
#define CGAL_TYPEDEF_H_

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_with_holes_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Epeck;
typedef CGAL::Polygon_2<Epeck> Polygon;
typedef CGAL::Polygon_with_holes_2<Epeck> Polygon_with_holes;
typedef CGAL::Aff_transformation_2<Epeck> Transformation;
typedef CGAL::Point_2<Epeck> Point;

#endif
