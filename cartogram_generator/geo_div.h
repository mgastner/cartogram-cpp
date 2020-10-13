#ifndef GEO_DIV_H_
#define GEO_DIV_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_with_holes_2<K> PolygonWH;

class GeoDiv {
  private:
    std::string id;
    std::vector<PolygonWH> polygons_with_holes;
  public:
    GeoDiv(const std::string);
    int n_polygons_with_holes(void) const;
    std::vector<PolygonWH> get_polygons_with_holes(void) const;
    void push_back(const PolygonWH);
};

#endif
