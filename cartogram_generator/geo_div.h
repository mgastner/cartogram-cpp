#ifndef GEO_DIV_H_
#define GEO_DIV_H_

#include "cgal_typedef.h"
#include <string>
#include <vector>

class GeoDiv {
  private:
    std::string id;
    std::vector<Polygon_with_holes> polygons_with_holes;
    GeoDiv();
  public:
    explicit GeoDiv(const std::string);
    int n_polygons_with_holes() const;
    std::vector<Polygon_with_holes> get_polygons_with_holes() const;
    std::vector<Polygon_with_holes> *ref_to_polygons_with_holes();
    void push_back(const Polygon_with_holes);
};

#endif
