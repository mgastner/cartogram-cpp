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

int PLL::pos() {
  return pos_;
}

Polyline PLL::pll() {
  return pll_;
}

Point PLL::v1() {
  return v1_;
}

Point PLL::v2() {
  return v2_;
}

Point PLL::vl() {
  return vl_;
}

Point PLL::v2l() {
  return v2l_;
}

int PLL::gd() {
  return gd_;
}

int PLL::pgnwh() {
  return pgnwh_;
}

bool PLL::is_hole() {
  return is_hole_;
}

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

void PLL::set_is_hole(bool b) {
  is_hole_ = b;
}
