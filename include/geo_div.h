#ifndef GEO_DIV_H_
#define GEO_DIV_H_

#include "intersection.h"
#include <string>
#include <vector>

class GeoDiv {
  private:
    std::set<std::string> adjacent_geodivs_;
    std::string id_;
    std::vector<Polygon_with_holes> polygons_with_holes_;
    GeoDiv();

  public:
    explicit GeoDiv(const std::string);
    const std::set<std::string> adjacent_geodivs() const;
    void adjacent_to(const std::string);
    double area() const;
    Bbox bbox() const;
    const std::string id() const;
    const std::vector<Segment> intersections(unsigned int) const;
    const Polygon_with_holes largest_polygon_with_holes() const;
    unsigned int n_points() const;
    unsigned int n_polygons_with_holes() const;
    unsigned int n_rings() const;
    Point point_on_surface_of_geodiv() const;
    Point point_on_surface_of_polygon_with_holes(
      const Polygon_with_holes) const;
    const std::vector<Polygon_with_holes> polygons_with_holes() const;
    void push_back(const Polygon_with_holes);
    std::vector<Polygon_with_holes> *ref_to_polygons_with_holes();
};

#endif
