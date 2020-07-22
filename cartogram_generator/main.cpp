#include <mapbox/geojson.hpp>
#include <mapbox/geojson/rapidjson.hpp>
#include <mapbox/geometry.hpp>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace mapbox::geojson;

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
}
