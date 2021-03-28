#include "geo_div.h"

GeoDiv::GeoDiv(const std::string i) : id_(i)
{
  return;
}

const std::string GeoDiv::id() const
{
  return id_;
}

int GeoDiv::n_polygons_with_holes() const
{
  return polygons_with_holes_.size();
}

const std::vector<Polygon_with_holes> GeoDiv::polygons_with_holes() const
{
  return polygons_with_holes_;
}

std::vector<Polygon_with_holes> *GeoDiv::ref_to_polygons_with_holes()
{
  return &polygons_with_holes_;
}

void GeoDiv::push_back(const Polygon_with_holes pgn_wh)
{
  polygons_with_holes_.push_back(pgn_wh);
  return;
}

double GeoDiv::area() const
{
  double a = 0.0;
  for (auto pwh : polygons_with_holes()) {
    Polygon ext_ring = pwh.outer_boundary();
    a += ext_ring.area();
    for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
      Polygon hole = *hci;
      a += hole.area();
      std::cout << "hole.area() = " << hole.area() << std::endl;
    }
  }
  return a;
}

XY_Point GeoDiv::centroid_of_polygon(const Polygon polyg) const
{
  // Code for centroid from: https://graphics.stanford.edu/courses/cs368-04-
  // spring/manuals/CGAL_Tutorial.pdf (accessed on 2021-Mar-15).
  // Check if the polygon has at least three vertices.
  assert(polyg.size() >= 3);
  Polygon::Vertex_circulator start = polyg.vertices_circulator();
  Polygon::Vertex_circulator cur = start;
  Polygon::Vertex_circulator next = cur;
  ++next;
  CGAL::Vector_2<Epick> centre(0, 0);
  double a = 0.0, atot = 0.0;
  do {
    a = ((*cur).x()) * ((*next).y()) - ((*next).x()) * ((*cur).y());
    centre = centre + a * ((*cur - CGAL::ORIGIN) + (*next - CGAL::ORIGIN));
    atot = atot + a;
    cur = next;
    ++next;
  } while (cur != start);
  atot = 3 * atot;
  centre /= atot;
  return XY_Point(CGAL::ORIGIN + centre);
}

XY_Point GeoDiv::centroid_of_polygon_with_holes(const Polygon_with_holes pwh)
const
{
  // Idea from https://math.stackexchange.com/questions/623841/finding-
  // centroid-of-a-polygon-with-holes
  Polygon ext_ring = pwh.outer_boundary();
  double a_ext = ext_ring.area();
  XY_Point c_ext = centroid_of_polygon(ext_ring);
  CGAL::Vector_2<Epick> weighted_sum_of_centroids(c_ext.x, c_ext.y);
  weighted_sum_of_centroids *= a_ext;
  double total_area = a_ext;
  for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
    Polygon hole = *hci;
    double a_hole = hole.area();
    XY_Point c_hole = centroid_of_polygon(hole);
    CGAL::Vector_2<Epick> c_hole_vector(c_hole.x, c_hole.y);
    weighted_sum_of_centroids += a_hole * c_hole_vector;
    total_area += a_hole;
  }
  weighted_sum_of_centroids /= total_area;
  XY_Point c_pwh(CGAL::ORIGIN + weighted_sum_of_centroids);
  std::cout << "Centroid of pwh: (" << c_pwh.x << ", " << c_pwh.y << ")" << std::endl;
  return c_pwh;
}

XY_Point GeoDiv::point_in_largest_polygon_with_holes() const
{
  // Find largest polygon with hole in GeoDiv
  double max_pwh_area = 0.0;
  unsigned int max_pwh_index = 0;
  for (int j = 0; j < n_polygons_with_holes(); ++j) {
    Polygon_with_holes pwh = polygons_with_holes()[j];
    Polygon ext_ring = pwh.outer_boundary();
    double a = ext_ring.area();
    for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
      Polygon hole = *hci;
      a += hole.area();
    }
    if (a > max_pwh_area) {
      max_pwh_area = a;
      max_pwh_index = j;
    }
  }
  Polygon_with_holes max_pwh = polygons_with_holes()[max_pwh_index];
  XY_Point centroid = centroid_of_polygon_with_holes(max_pwh);
  if (centroid.is_in_polygon_with_holes(max_pwh)) {
    std::cout << "Nominally in pwh" << std::endl;
  }
  return centroid;
}
