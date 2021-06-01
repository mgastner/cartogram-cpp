#include <string>
#include <vector>

#include "pll.h"

PLL::PLL(int pos,
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
int PLL::pos() {
  return pos_;
}

// Get the CGAL polyline object stored inside
Polyline PLL::pll() {
  return pll_;
}

// Get its first vertex
Point PLL::v1() {
  return v1_;
}

// Get its second vertex
Point PLL::v2() {
  return v2_;
}

// Get its last vertex
Point PLL::vl() {
  return vl_;
}

// Get its second last vertex
Point PLL::v2l() {
  return v2l_;
}

// Get its associated GeoDiv
int PLL::gd() {
  return gd_;
}

// Get its associated Polygon_with_holes
int PLL::pgnwh() {
  return pgnwh_;
}

// Get whether it is a hole or not
bool PLL::is_hole() {
  return is_hole_;
}

// Set/assign the CGAL polyline object stored inside
void PLL::set_pll(Polyline pll) {
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
void PLL::set_is_hole(bool b) {
  is_hole_ = b;
}
