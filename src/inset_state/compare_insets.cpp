#include "compare_insets.h"
#include <list>
#include <CGAL/centroid.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <boost/geometry.hpp>
#include <boost/geometry/algorithms/discrete_frechet_distance.hpp>
#include <boost/geometry/algorithms/discrete_hausdorff_distance.hpp>
#include <boost/geometry/algorithms/sym_difference.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/multi_linestring.hpp>

// https://stackoverflow.com/questions/38810048/intersection-of-boostgeometrymodellinestring-with-boostgeometrymodel
using point_type = boost::geometry::model::d2::point_xy<double>;
using linestring_type = boost::geometry::model::linestring<point_type>;
using multi_linestring_type =
  boost::geometry::model::multi_linestring<linestring_type>;

bool account_for_holes = false;

// Debugging test for checking the order of line strings in the Hausdorff
// distance function. The first line string in the argument is the anchor
// a.k.a. the line string of which we take each point and measure its distance
// to the closest point in the second line string.
// void hausdorff_test() {

//   linestring_type ls1, ls2;
//   boost::geometry::read_wkt("LINESTRING(1 1,2 1,3 1,4 1,5 1)", ls1);
//   boost::geometry::read_wkt("LINESTRING(1 2,2 4,3 2,4 2,5 2)", ls2);

//   double res1 = boost::geometry::discrete_hausdorff_distance(ls1, ls2);
//   double res2 = boost::geometry::discrete_hausdorff_distance(ls2, ls1);

//   std::cout << "Discrete Hausdorff Distance 1: " << res1 << std::endl;
//   std::cout << "Discrete Hausdorff Distance 2: " << res2 << std::endl;
//   return;
// }

double area_polygon_with_holes(Polygon_with_holes poly) {
  double a = 0;
  a += poly.outer_boundary().area();
  for (auto h = poly.holes_begin(); h != poly.holes_end(); ++h) {
    a += h->area();
  }
  return a;
}

Polygon_with_holes normalised_polygon(Polygon_with_holes poly) {

  // How to calculate areas when normalising polygons
  double poly_area;
  if (account_for_holes){
    poly_area = area_polygon_with_holes(poly);
  } else {
    poly_area = poly.outer_boundary().area();
  }

  // Get the scale factor and initialise a CGAL scale object.
  double scale_factor = sqrt(1 / poly_area);
  const Transformation scale(CGAL::SCALING, scale_factor);

  // Construct the new polygon.
  Polygon new_outer =
    transform(scale, poly.outer_boundary());
  std::vector<Polygon> new_holes;
  for (auto h = poly.holes_begin(); h != poly.holes_end(); ++h) {
    new_holes.push_back(transform(scale, *h));
  }
  Polygon_with_holes new_poly(new_outer, new_holes.begin(), new_holes.end());
  return new_poly;
}

Polygon_with_holes translated_polygon (Polygon_with_holes poly,
                                       double x_shift,
                                       double y_shift) {

  // Get the transformation.
  const Transformation translate(CGAL::TRANSLATION,
                                 CGAL::Vector_2<Scd> (x_shift, y_shift));

  // Construct the new polygon.
  Polygon new_outer =
    transform(translate, poly.outer_boundary());
  std::vector<Polygon> new_holes;
  for (auto h = poly.holes_begin(); h != poly.holes_end(); ++h) {
    new_holes.push_back(transform(translate, *h));
  }
  Polygon_with_holes new_poly(new_outer, new_holes.begin(), new_holes.end());
  return new_poly;
}

std::pair<double, double> get_translation_coordinates (
  Polygon_with_holes anchor_poly,
  Polygon_with_holes displaced_poly)
{

  // Note: The centroid calculation ignores the holes and returns
  // a polygon's cenroid purely on the basis of its outer boundary.
  Point anchor_point = CGAL::centroid(
    anchor_poly.outer_boundary().begin(),
    anchor_poly.outer_boundary().end(),
    CGAL::Dimension_tag<0>()
  );
  Point displaced_point = CGAL::centroid(
    displaced_poly.outer_boundary().begin(),
    displaced_poly.outer_boundary().end(),
    CGAL::Dimension_tag<0>()
  );

  std::pair<double, double> displacement;
  displacement.first = anchor_point.x() - displaced_point.x();
  displacement.second = anchor_point.y() - displaced_point.y();

  return displacement;
}

linestring_type get_linestring(Polygon poly) {
  linestring_type orig_outer_linestring;
  for (auto pt : poly) {
    boost::geometry::append(
      orig_outer_linestring, point_type(pt.x(), pt.y())
    );
  }
  return orig_outer_linestring;
}

// Debug function to test out Frechet distance variations.
Polygon_with_holes get_polygon_with_new_starting_point (
  unsigned int start_pt,
  Polygon_with_holes poly)
{
  std::vector<Point> start_points;
  Polygon outer_pwh = poly.outer_boundary();
  Polygon new_outer_pwh;
  for (unsigned int i = 0; i < outer_pwh.size(); ++i) {
    if (i < start_pt) {
      start_points.push_back(outer_pwh[i]);
    } else {
      new_outer_pwh.push_back(outer_pwh[i]);
    }
  }
  for (auto pt : start_points) {
    new_outer_pwh.push_back(pt);
  }
  return Polygon_with_holes (
    new_outer_pwh, poly.holes_begin(), poly.holes_end()
  );
}

double frechet_distance(Polygon_with_holes orig_poly,
                        Polygon_with_holes new_poly) {
  double distance = 0;
  linestring_type orig_outer_linestring =
    get_linestring(orig_poly.outer_boundary());
  linestring_type new_outer_linestring =
    get_linestring(new_poly.outer_boundary());
  distance = std::max(
    distance,
    boost::geometry::discrete_frechet_distance(
      orig_outer_linestring, new_outer_linestring
    )
  );

  if (account_for_holes){
    std::deque<Polygon> orig_holes = orig_poly.holes();
    std::deque<Polygon> new_holes = new_poly.holes();
    for(unsigned int i = 0; i < orig_holes.size(); ++i) {
      linestring_type orig_hole_linestring = get_linestring(orig_holes[i]);
      linestring_type new_hole_linestring = get_linestring(new_holes[i]);
      distance = std::max(
        distance,
        boost::geometry::discrete_frechet_distance(
          orig_hole_linestring, new_hole_linestring
        )
      );
    }
  }
  return distance;
}

double hausdorff_distance(Polygon_with_holes orig_poly,
                          Polygon_with_holes new_poly) {
  double distance = 0;
  linestring_type orig_outer_linestring =
    get_linestring(orig_poly.outer_boundary());
  linestring_type new_outer_linestring =
    get_linestring(new_poly.outer_boundary());
  distance = std::max(
    distance,
    boost::geometry::discrete_hausdorff_distance(
      orig_outer_linestring, new_outer_linestring
    )
  );

  if (account_for_holes){
    std::deque<Polygon> orig_holes = orig_poly.holes();
    std::deque<Polygon> new_holes = new_poly.holes();
    for(unsigned int i = 0; i < orig_holes.size(); ++i) {
      linestring_type orig_hole_linestring = get_linestring(orig_holes[i]);
      linestring_type new_hole_linestring = get_linestring(new_holes[i]);
      distance = std::max(
        distance,
        boost::geometry::discrete_hausdorff_distance(
          orig_hole_linestring, new_hole_linestring
        )
      );
    }
  }
  return distance;
}

double symmetric_distance(Polygon_with_holes orig_poly,
                          Polygon_with_holes new_poly) {

  // Calculate the intersection of the original and transformed polygons.
  std::list<Polygon_with_holes> inter_polys;
  if (account_for_holes) {
    CGAL::intersection(orig_poly, new_poly, std::back_inserter(inter_polys));
  } else {
    CGAL::intersection(
      orig_poly.outer_boundary(),
      new_poly.outer_boundary(),
      std::back_inserter(inter_polys)
    );
  }

  double inter_area = 0;
  for (auto poly : inter_polys){
    inter_area += area_polygon_with_holes(poly);
  }
  double sum_orig_new_area;
  if (account_for_holes) {
    sum_orig_new_area =
      area_polygon_with_holes(orig_poly) + area_polygon_with_holes(new_poly);
  } else {
    sum_orig_new_area =
      orig_poly.outer_boundary().area() + new_poly.outer_boundary().area();
  }
  return (sum_orig_new_area - inter_area * 2) /
    (sum_orig_new_area - inter_area);
}

std::map<std::string, double> InsetState::get_geo_div_differences(
  std::function<double (Polygon_with_holes, Polygon_with_holes)> distance_func,
  bool compare_with_original)
{
  
  // Initialize a vector to store the max distance for each GeoDiv
  std::map<std::string, double> differences_geo_divs;

  const std::vector<GeoDiv> *original_geo_divs;
  const std::vector<GeoDiv> *new_geo_divs;
  if (compare_with_original) {
    original_geo_divs = &geo_divs_original_;
    new_geo_divs = &geo_divs_;
  } else {
    original_geo_divs = &geo_divs_;
    new_geo_divs = &geo_divs_compare_;
  }

  // Iterate through each GeoDiv
  for (unsigned int i = 0; i < original_geo_divs -> size(); ++i) {

    // Vectors of the current GeoDiv's polygon with holes,
    // before and after integration
    std::vector<Polygon_with_holes> original_pwh_vector =
      (*original_geo_divs)[i].polygons_with_holes();
    std::vector<Polygon_with_holes> new_pwh_vector =
      (*new_geo_divs)[i].polygons_with_holes();

    // Iterate through each pair of polygon with holes.
    std::vector<double> distances;
    for (unsigned int j = 0; j < original_pwh_vector.size(); ++j) {

      // Normalise polygons (outer boundary area = 1)
      Polygon_with_holes norm_original_pwh =
        normalised_polygon(original_pwh_vector[j]);
      Polygon_with_holes norm_new_pwh =
        normalised_polygon(new_pwh_vector[j]);

      // Move the new polygon to have its centroid coincide with that of
      // the original (pre-integration) polygon
      std::pair<double, double> translation =
        get_translation_coordinates(norm_original_pwh, norm_new_pwh);
      Polygon_with_holes new_pwh =
        translated_polygon(
          norm_new_pwh, translation.first, translation.second
        );
      
      // DEBUG: Change outer boundary starting point of one of the polygons
      // to check Frechet distance change.
      // new_pwh = get_polygon_with_new_starting_point(20, new_pwh);
      
      // Calculate distance
      distances.push_back(distance_func(norm_original_pwh, new_pwh));
    }

    // Get the maximum distance and take that as the "distance" of the whole
    // GeoDiv
    double dist = *std::max_element(distances.begin(), distances.end());
    std::string geo_div_id = (*new_geo_divs)[i].id();
    differences_geo_divs.insert(
      std::pair<std::string, double> (geo_div_id, dist)
    );
  }
  return differences_geo_divs;
}
