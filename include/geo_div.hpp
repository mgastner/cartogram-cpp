#ifndef GEO_DIV_HPP_
#define GEO_DIV_HPP_

#include "ellipse.hpp"
#include "pwh.hpp"

class GeoDiv
{
private:
  std::set<std::string> adjacent_geodivs_;
  std::string id_;
  std::vector<Ellipse> min_ellipses_;
  std::vector<Polygon_with_holes> polygons_with_holes_;
  GeoDiv();

public:
  explicit GeoDiv(std::string);
  [[nodiscard]] const std::set<std::string> &adjacent_geodivs() const;
  void adjacent_to(const std::string &);
  [[nodiscard]] double area() const;
  void clear_min_ellipses();
  [[nodiscard]] Bbox bbox() const;
  [[nodiscard]] const std::string &id() const;
  [[nodiscard]] Polygon_with_holes largest_polygon_with_holes() const;
  const std::vector<Ellipse> &min_ellipses() const;
  [[nodiscard]] size_t n_points() const;
  [[nodiscard]] size_t n_polygons_with_holes() const;
  [[nodiscard]] unsigned int n_rings() const;
  [[nodiscard]] Point point_on_surface_of_geodiv() const;
  [[nodiscard]] Point point_on_surface_of_polygon_with_holes(
    const Polygon_with_holes &) const;
  [[nodiscard]] const std::vector<Polygon_with_holes> &polygons_with_holes()
    const;
  void push_back(const Ellipse &);
  void push_back(const Polygon_with_holes &);
  std::vector<Polygon_with_holes> &ref_to_polygons_with_holes();
  void sort_pwh_descending_by_area();
  void update_id(const std::string &);
};

#endif  // GEO_DIV_HPP_
