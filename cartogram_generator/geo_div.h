#ifndef GEO_DIV_H_
#define GEO_DIV_H_

#include "cgal_typedef.h"
#include <string>
#include <vector>

class GeoDiv {
  private:
    std::string id;
    std::vector<Polygon_with_holes> polygons_with_holes;
  public:
    void set_id(const std::string);
    int n_polygons_with_holes(void) const;
    std::vector<Polygon_with_holes> get_polygons_with_holes(void) const;
    std::vector<Polygon_with_holes> *ref_to_polygons_with_holes(void);
    void push_back(const Polygon_with_holes);
};

#endif
