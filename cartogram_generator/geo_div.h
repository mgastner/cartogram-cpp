#ifndef GEO_DIV_H_
#define GEO_DIV_H_

#include "cgal_typedef.h"
#include "xy_point.h"
#include <string>
#include <vector>

class GeoDiv {
  private:
    std::string id_;
    std::vector<Polygon_with_holes> polygons_with_holes_;
    GeoDiv();
  public:
    explicit GeoDiv(const std::string);
    const std::string id() const;
    int n_polygons_with_holes() const;
    const std::vector<Polygon_with_holes> polygons_with_holes() const;
    std::vector<Polygon_with_holes> *ref_to_polygons_with_holes();
    void push_back(const Polygon_with_holes);
    double area() const;
    XY_Point centroid_of_polygon(const Polygon) const;
    XY_Point centroid_of_polygon_with_holes(const Polygon_with_holes) const;
    XY_Point point_on_surface(const Polygon_with_holes) const;
    XY_Point point_in_largest_polygon_with_holes() const;
};

#endif
