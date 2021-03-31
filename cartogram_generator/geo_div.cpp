#include "constants.h"
#include "geo_div.h"
#include "find_intersections_at_y.h"

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

// Function that takes a Polygon_with_holes and returns the midpoint of the
// longest line segment that is inside the polygon and is halfway through the
// northern and southern tip of the polygon. Reference:
// https://gis.stackexchange.com/questions/76498/how-is-st-pointonsurface-
// calculated
XY_Point GeoDiv::point_on_surface(const Polygon_with_holes pwh) const
{
  // Calculate line_y
  CGAL::Bbox_2 bb = pwh.bbox();
  double line_y = (bb.ymin() + bb.ymax()) / 2;

  // Epsilon based on default resolution
  double epsilon = 1e-6 * (1.0/sub_sample_resolution);

  // Vector to store intersections
  std::vector<intersection> intersections;

  // Getting outer_boundary from pwh
  Polygon ext_ring = pwh.outer_boundary();

  // Setting up previous point to form segment (curr_point, prev_point)
  XY_Point prev_point;
  prev_point.x = ext_ring[ext_ring.size()-1][0];
  prev_point.y = ext_ring[ext_ring.size()-1][1];

  // Finding all the intersections with exterior ring
  for (unsigned int l = 0; l < ext_ring.size(); ++l) {
    XY_Point curr_point;
    curr_point.x = ext_ring[l][0];
    curr_point.y = ext_ring[l][1];
    intersection temp;

    // Function to calculate whether intersection exists between
    // segment (prev_point, curr_point) and line_y
    if (line_y_intersects(curr_point,
                          prev_point,
                          line_y,
                          &temp,
                          0,  // no target_density required
                          epsilon)) {
      intersections.push_back(temp);
    }
    prev_point.x = curr_point.x;
    prev_point.y = curr_point.y;
  }

  // Finding all intersections for holes
  for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
    Polygon hole = *hci;
    prev_point.x = hole[hole.size()-1][0];
    prev_point.y = hole[hole.size()-1][1];
    for (unsigned int l = 0; l < hole.size(); ++l) {
      XY_Point curr_point;
      curr_point.x = hole[l][0];
      curr_point.y = hole[l][1];
      intersection temp;
      if (line_y_intersects(curr_point,
                            prev_point,
                            line_y,
                            &temp,
                            0,
                            epsilon)) {
        intersections.push_back(temp);
      }
      prev_point.x = curr_point.x;
      prev_point.y = curr_point.y;
    }
  }

  // Check if the number of intersections is odd
  if (intersections.size() % 2 != 0) {
    std::cerr << "Unable to calculate point on surface!" << std::endl;
    std::cerr << "Incorrect Topology" << std::endl;
    std::cerr << "Number of intersections: " << intersections.size();
    std::cerr << std::endl;
    std::cerr << "Y-coordinate: " << line_y << std::endl;
    std::cerr << "Intersection points: " << std::endl;
    for (unsigned int l = 0; l < intersections.size(); ++l) {
      std::cerr << intersections[l].x << std::endl;
    }
    _Exit(932874);
  }
  std::sort(intersections.begin(), intersections.end());

  // Assign directions (i.e. whether the line is entering or leaving the
  // polygon with holes)
  for (unsigned int l = 0; l < intersections.size(); ++l) {
    intersections[l].direction = (l%2 == 0);
  }

  // Assigning length of line segments (to find longest) using the
  // target_density property of intersections for line segment lengths
  for (unsigned int l = 0; l < intersections.size(); l += 2) {
    intersections[l].target_density =
      intersections[l + 1].x - intersections[l].x;
  }

  // Finding maximum segment length
  double max_length = intersections[0].target_density;
  double left = intersections[0].x;
  double right = intersections[1].x;
  XY_Point midpoint((right + left) / 2, line_y);

  // Iterating through lengths
  for (unsigned int l = 0; l < intersections.size(); l += 2) { \
    if (intersections[l].target_density >= max_length) {
      left = intersections[l].x;
      right = intersections[l + 1].x;
      max_length = intersections[l].target_density;
      midpoint.x = (right + left) / 2;
    }
  }

  // Making final midpoint and returning it
  return midpoint;
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
  return centroid.is_in_polygon_with_holes(max_pwh) ? centroid :
         point_on_surface(max_pwh);
}
