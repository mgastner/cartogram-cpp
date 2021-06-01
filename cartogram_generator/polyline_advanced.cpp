#include <string>
#include <vector>

#include "polyline_advanced.h"

Polyline_advanced::Polyline_advanced(int pos,
    Polyline pll,
    int gd,
    int pgnwh,
    bool is_hole) {
  pos_ = pos;
  pll_ = pll;
  gd_ = gd;
  pgnwh_ = pgnwh;

  v1_ = *pll.begin();
  auto v2_it = pll.begin();
  v2_it++;
  v2_ = *v2_it;
  vl_ = *pll.rbegin();
  auto v2l_it = pll.rbegin();
  v2l_it++;
  v2l_ = *v2l_it;

  is_hole_ = is_hole;
}

// Get the position of the polyline inside ct_polylines
int Polyline_advanced::pos() {
  return pos_;
}

// Get the CGAL polyline object stored inside
Polyline Polyline_advanced::pll() {
  return pll_;
}

// Get its first vertex
Point Polyline_advanced::v1() {
  return v1_;
}

// Get its second vertex
Point Polyline_advanced::v2() {
  return v2_;
}

// Get its last vertex
Point Polyline_advanced::vl() {
  return vl_;
}

// Get its second last vertex
Point Polyline_advanced::v2l() {
  return v2l_;
}

// Get its associated GeoDiv
int Polyline_advanced::gd() {
  return gd_;
}

// Get its associated Polygon_with_holes
int Polyline_advanced::pgnwh() {
  return pgnwh_;
}

// Get whether it is a hole or not
bool Polyline_advanced::is_hole() {
  return is_hole_;
}

// Set/assign the CGAL polyline object stored inside
void Polyline_advanced::set_pll(Polyline pll) {
  pll_ = pll;
  v1_ = *pll.begin();
  auto v2_it = pll.begin();
  v2_it++;
  v2_ = *v2_it;
  vl_ = *pll.rbegin();
  auto v2l_it = pll.rbegin();
  v2l_it++;
  v2l_ = *v2l_it;
}

// Set/assign whether it is a hole or not
void Polyline_advanced::set_is_hole(bool b) {
  is_hole_ = b;
}
