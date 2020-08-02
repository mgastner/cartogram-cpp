#include <mapbox/geojson.hpp>
#include <mapbox/geojson/rapidjson.hpp>
#include <mapbox/geometry.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace mapbox::geojson;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;

template <typename T = geojson>
geojson readGeoJSON(const string &path, bool use_convert) {
  ifstream t(path.c_str());
  stringstream buffer;
  buffer << t.rdbuf();
  if (use_convert) {
    rapidjson_document d;
    d.Parse<0>(buffer.str().c_str());
    return convert<T>(d);
  } else {
    return parse(buffer.str());
  }
}
int main() {
  bool use_convert = true;
  const auto data = readGeoJSON("../sample_data/feature-collection.json",
                                use_convert);
  //const auto &features = data.get<feature_collection>();
  Point points[] = { Point(0,0), Point(1,0), Point(0,1), Point(1,1)};
  Polygon_2 pgn(points, points+4);
  // check if the polygon is simple.
  cout << "The polygon is " <<
    (pgn.is_simple() ? "" : "not ") << "simple." << endl;
  // check if the polygon is convex
  cout << "The polygon is " <<
    (pgn.is_convex() ? "" : "not ") << "convex." << endl;

  return 0;
}
